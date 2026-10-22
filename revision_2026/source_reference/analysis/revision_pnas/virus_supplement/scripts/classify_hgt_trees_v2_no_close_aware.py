#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automatic phylogenetic classification for large-scale HGT candidate trees.

Design goals
------------
1. Do not force a binary HGT / non-HGT decision.
2. Keep most results, but assign evidence levels and flags.
3. Treat close-rich trees differently from close-sparse trees.
4. Mark possible long-branch, ORF/frame, alignment, and donor/self mixed topology issues.
5. Use only Python standard library.

Expected input tree labels
--------------------------
Leaf labels are expected to start with one of:
  QUERY|..., SELF|..., CLOSE|..., DONOR_OUT|..., OUT_..., UNKNOWN|...
The script is intentionally tolerant and will infer categories from prefixes.

IQ-TREE support labels
----------------------
Internal node labels like 98.6/100 are interpreted as:
  SH-aLRT / UFBoot
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set, Iterable


# -----------------------------
# Basic Newick parser
# -----------------------------

class Node:
    __slots__ = ("name", "length", "children", "parent", "sh_alrt", "ufboot", "leafset")

    def __init__(self, name: str = "", length: float = 0.0):
        self.name = name
        self.length = length
        self.children: List[Node] = []
        self.parent: Optional[Node] = None
        self.sh_alrt: Optional[float] = None
        self.ufboot: Optional[float] = None
        self.leafset: Optional[Set[str]] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children


def _parse_support(label: str) -> Tuple[Optional[float], Optional[float]]:
    label = (label or "").strip()
    if not label:
        return None, None
    # IQ-TREE often writes SH-aLRT/UFBoot as "98.6/100".
    m = re.match(r"^\s*([0-9.]+)\s*/\s*([0-9.]+)\s*$", label)
    if m:
        try:
            return float(m.group(1)), float(m.group(2))
        except Exception:
            return None, None
    # Sometimes a single bootstrap value is present; treat as UFBoot only.
    m = re.match(r"^\s*([0-9.]+)\s*$", label)
    if m:
        try:
            return None, float(m.group(1))
        except Exception:
            return None, None
    return None, None


def parse_newick(text: str) -> Node:
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    n = len(s)
    i = 0

    def skip_ws(j: int) -> int:
        while j < n and s[j].isspace():
            j += 1
        return j

    def parse_label(j: int) -> Tuple[str, int]:
        j = skip_ws(j)
        start = j
        while j < n and s[j] not in ":,()":
            j += 1
        return s[start:j].strip(), j

    def parse_length(j: int) -> Tuple[float, int]:
        j = skip_ws(j)
        if j >= n or s[j] != ":":
            return 0.0, j
        j += 1
        start = j
        while j < n and s[j] not in ",()":
            j += 1
        raw = s[start:j].strip()
        try:
            return float(raw), j
        except Exception:
            return 0.0, j

    def parse_subtree(j: int) -> Tuple[Node, int]:
        j = skip_ws(j)
        if j < n and s[j] == "(":
            j += 1
            node = Node()
            while True:
                child, j = parse_subtree(j)
                child.parent = node
                node.children.append(child)
                j = skip_ws(j)
                if j < n and s[j] == ",":
                    j += 1
                    continue
                if j < n and s[j] == ")":
                    j += 1
                    break
                raise ValueError(f"Malformed Newick near position {j}: {s[j:j+60]!r}")
            label, j = parse_label(j)
            node.name = label
            node.sh_alrt, node.ufboot = _parse_support(label)
            node.length, j = parse_length(j)
            return node, j
        else:
            label, j = parse_label(j)
            node = Node(label)
            node.length, j = parse_length(j)
            return node, j

    root, i = parse_subtree(0)
    i = skip_ws(i)
    if i != n:
        # Be tolerant but warn through exception for clearly malformed input.
        rest = s[i:i+80]
        if rest.strip():
            raise ValueError(f"Unexpected trailing Newick text at {i}: {rest!r}")
    root.parent = None
    return root


def iter_nodes(root: Node) -> Iterable[Node]:
    stack = [root]
    while stack:
        x = stack.pop()
        yield x
        stack.extend(reversed(x.children))


def leaves(root: Node) -> List[Node]:
    return [x for x in iter_nodes(root) if x.is_leaf]


def postorder(root: Node) -> Iterable[Node]:
    for c in root.children:
        yield from postorder(c)
    yield root


# -----------------------------
# Category / metadata utilities
# -----------------------------

def leaf_category(label: str) -> str:
    x = (label or "").strip()
    prefix = x.split("|", 1)[0].upper()
    if prefix == "QUERY" or x.startswith("QUERY"):
        return "query"
    if prefix == "SELF":
        return "self"
    if prefix == "CLOSE":
        return "close"
    # Treat all distant/out donor classes as donor-like.
    if prefix.startswith("DONOR") or prefix.startswith("OUT") or prefix.startswith("DISTANT"):
        return "donor"
    if prefix in {"UNKNOWN", "UNK"}:
        return "unknown"
    # Fallback: infer from substring.
    ux = x.upper()
    if ux.startswith("SELF|"):
        return "self"
    if ux.startswith("CLOSE|"):
        return "close"
    if ux.startswith("DONOR") or "|DONOR" in ux or ux.startswith("OUT_"):
        return "donor"
    return "other"


def parse_hgt_id_from_treefile(path: Path) -> str:
    # job dir name is the most reliable HGT id.
    return path.parent.name


def read_tsv_by_key(path: Optional[Path], key_col: str = "hgt_id") -> Dict[str, Dict[str, str]]:
    if not path or not path.exists():
        return {}
    out: Dict[str, Dict[str, str]] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or key_col not in reader.fieldnames:
            return {}
        for row in reader:
            k = row.get(key_col, "")
            if k:
                out[k] = {kk: (vv if vv is not None else "") for kk, vv in row.items()}
    return out


def safe_float(x: object, default: float = math.inf) -> float:
    try:
        if x is None or x == "":
            return default
        return float(x)
    except Exception:
        return default


def safe_int(x: object, default: int = 0) -> int:
    try:
        if x is None or x == "":
            return default
        return int(float(str(x)))
    except Exception:
        return default


# -----------------------------
# Tree metrics
# -----------------------------

def build_adjacency(root: Node) -> Dict[int, List[Tuple[int, float]]]:
    adj: Dict[int, List[Tuple[int, float]]] = defaultdict(list)
    for node in iter_nodes(root):
        if node.parent is not None:
            a, b = id(node), id(node.parent)
            w = max(0.0, float(node.length or 0.0))
            adj[a].append((b, w))
            adj[b].append((a, w))
    return adj


def dijkstra_from_query(query: Node, root: Node) -> Dict[int, float]:
    adj = build_adjacency(root)
    src = id(query)
    dist = {src: 0.0}
    heap = [(0.0, src)]
    while heap:
        d, u = heapq.heappop(heap)
        if d != dist.get(u):
            continue
        for v, w in adj.get(u, []):
            nd = d + w
            if nd < dist.get(v, math.inf):
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return dist


def support_level(sh: Optional[float], uf: Optional[float], sh_strong: float, uf_strong: float,
                  sh_mod: float, uf_mod: float) -> str:
    has = (sh is not None) or (uf is not None)
    if not has:
        return "none"
    sh_ok_strong = True if sh is None else sh >= sh_strong
    uf_ok_strong = True if uf is None else uf >= uf_strong
    if sh_ok_strong and uf_ok_strong:
        return "strong"
    sh_ok_mod = True if sh is None else sh >= sh_mod
    uf_ok_mod = True if uf is None else uf >= uf_mod
    if sh_ok_mod and uf_ok_mod:
        return "moderate"
    return "low"


def collect_leafsets(root: Node) -> None:
    for node in postorder(root):
        if node.is_leaf:
            node.leafset = {node.name}
        else:
            st: Set[str] = set()
            for c in node.children:
                st.update(c.leafset or set())
            node.leafset = st


def find_best_query_split(root: Node, query_name: str, max_local_split: int) -> Dict[str, object]:
    """Find the smallest internal edge split whose query-side has query plus at least one leaf.

    This is more root-tolerant than simply using query.parent, because if the tree is
    arbitrarily rooted next to QUERY, the complement side of some internal node may
    still represent the local query grouping.
    """
    collect_leafsets(root)
    all_leaves = root.leafset or set()
    best = None
    for node in iter_nodes(root):
        if node.is_leaf or node is root:
            continue
        clade = node.leafset or set()
        if not clade or clade == all_leaves:
            continue
        q_side = clade if query_name in clade else (all_leaves - clade)
        if query_name not in q_side:
            continue
        size = len(q_side)
        if size < 2 or size > max_local_split:
            continue
        cats = Counter(leaf_category(x) for x in q_side if x != query_name)
        sh, uf = node.sh_alrt, node.ufboot
        rec = {
            "query_split_size": size,
            "query_split_self": cats.get("self", 0),
            "query_split_close": cats.get("close", 0),
            "query_split_donor": cats.get("donor", 0),
            "query_split_unknown": cats.get("unknown", 0) + cats.get("other", 0),
            "query_split_sh_alrt": sh,
            "query_split_ufboot": uf,
            "query_split_support_label": node.name,
        }
        # Prefer smallest query-side; tie-break by support availability/level.
        if best is None:
            best = rec
        else:
            if size < best["query_split_size"]:
                best = rec
            elif size == best["query_split_size"]:
                old_sup = (best.get("query_split_sh_alrt") is not None) + (best.get("query_split_ufboot") is not None)
                new_sup = (sh is not None) + (uf is not None)
                if new_sup > old_sup:
                    best = rec
    if best is None:
        return {
            "query_split_size": "NA",
            "query_split_self": 0,
            "query_split_close": 0,
            "query_split_donor": 0,
            "query_split_unknown": 0,
            "query_split_sh_alrt": "NA",
            "query_split_ufboot": "NA",
            "query_split_support_label": "NA",
        }
    return best


def summarize_tree(treefile: Path, args, orf_meta: Dict[str, str]) -> Dict[str, object]:
    hgt_id = parse_hgt_id_from_treefile(treefile)
    text = treefile.read_text().strip()
    root = parse_newick(text)
    lf = leaves(root)
    if not lf:
        raise ValueError("no leaves")

    cats = {x.name: leaf_category(x.name) for x in lf}
    query_leaves = [x for x in lf if cats[x.name] == "query"]
    if not query_leaves:
        raise ValueError("no QUERY leaf")
    query = query_leaves[0]
    query_name = query.name

    counts = Counter(cats.values())
    counts.pop("query", None)
    recipient_count = counts.get("self", 0) + counts.get("close", 0)

    terminal_lengths = [max(0.0, float(x.length or 0.0)) for x in lf if x is not query]
    median_terminal = statistics.median(terminal_lengths) if terminal_lengths else 0.0
    query_branch = max(0.0, float(query.length or 0.0))
    query_branch_ratio = query_branch / median_terminal if median_terminal > 0 else math.inf

    dist = dijkstra_from_query(query, root)
    leaf_dist = []
    min_dist_by_cat = defaultdict(lambda: math.inf)
    for x in lf:
        if x is query:
            continue
        cat = cats[x.name]
        d = dist.get(id(x), math.inf)
        leaf_dist.append((d, cat, x.name))
        if d < min_dist_by_cat[cat]:
            min_dist_by_cat[cat] = d
    leaf_dist.sort(key=lambda z: z[0])

    nearest_k = max(1, args.nearest_k)
    near = leaf_dist[:nearest_k]
    near5 = leaf_dist[:5]
    near10 = leaf_dist[:10]
    near_counts = Counter(cat for _, cat, _ in near)
    near5_counts = Counter(cat for _, cat, _ in near5)
    near10_counts = Counter(cat for _, cat, _ in near10)

    d_self = min_dist_by_cat.get("self", math.inf)
    d_close = min_dist_by_cat.get("close", math.inf)
    d_recipient = min(d_self, d_close)
    d_donor = min_dist_by_cat.get("donor", math.inf)
    d_unknown = min(min_dist_by_cat.get("unknown", math.inf), min_dist_by_cat.get("other", math.inf))

    # Root-tolerant local split around QUERY.
    split = find_best_query_split(root, query_name, args.max_local_split)
    sh = split.get("query_split_sh_alrt")
    uf = split.get("query_split_ufboot")
    sh_f = sh if isinstance(sh, float) else None
    uf_f = uf if isinstance(uf, float) else None
    sup_level = support_level(sh_f, uf_f, args.sh_strong, args.uf_strong, args.sh_moderate, args.uf_moderate)

    # ORF/frame metadata.
    stop_codons = safe_int(orf_meta.get("stop_codons_to_X", orf_meta.get("stop_codons", "")), default=-1)
    amb_codons = safe_int(orf_meta.get("ambiguous_codons_to_X", orf_meta.get("ambiguous_codons", "")), default=-1)
    orf_status = orf_meta.get("orf_status", orf_meta.get("reason", ""))
    frame_status = orf_meta.get("status", "")
    best_strand = orf_meta.get("best_strand", orf_meta.get("selected_strand", ""))
    best_frame = orf_meta.get("best_frame_0based", orf_meta.get("selected_frame_0based", ""))

    orf_known = bool(orf_meta)
    orf_bad = False
    orf_minor = False
    orf_clean = False
    if orf_known:
        if stop_codons >= args.bad_stop_cutoff or "internal" in orf_status.lower():
            orf_bad = True
        elif stop_codons >= args.minor_stop_cutoff or "minor" in orf_status.lower():
            orf_minor = True
        if stop_codons == 0 or "clean" in orf_status.lower():
            orf_clean = True

    long_branch = ((query_branch >= args.long_branch_abs and query_branch_ratio >= args.long_branch_ratio)
                   or query_branch >= args.very_long_branch_abs)

    donor_nearer = math.isfinite(d_donor) and math.isfinite(d_recipient) and (d_donor < d_recipient * args.distance_ratio)
    recipient_nearer = math.isfinite(d_donor) and math.isfinite(d_recipient) and (d_recipient < d_donor * args.distance_ratio)
    donor_only_context = math.isfinite(d_donor) and not math.isfinite(d_recipient)
    recipient_only_context = math.isfinite(d_recipient) and not math.isfinite(d_donor)

    nearest1_cat = leaf_dist[0][1] if leaf_dist else "NA"
    nearest1_name = leaf_dist[0][2] if leaf_dist else "NA"
    nearest1_dist = leaf_dist[0][0] if leaf_dist else math.inf

    near10_recipient = near10_counts.get("self", 0) + near10_counts.get("close", 0)
    near10_donor = near10_counts.get("donor", 0)
    near5_recipient = near5_counts.get("self", 0) + near5_counts.get("close", 0)
    near5_donor = near5_counts.get("donor", 0)

    split_recipient = safe_int(split.get("query_split_self"), 0) + safe_int(split.get("query_split_close"), 0)
    split_donor = safe_int(split.get("query_split_donor"), 0)
    local_mixed = split_recipient > 0 and split_donor > 0
    nearest_mixed = near10_recipient > 0 and near10_donor > 0

    close_count = counts.get("close", 0)
    self_count = counts.get("self", 0)
    donor_count = counts.get("donor", 0)
    close_rich = close_count >= args.min_close_rich
    close_sparse = close_count <= args.few_close
    recipient_rich = recipient_count >= args.min_recipient_rich

    flags: List[str] = []
    if close_count == 0:
        flags.append("no_close_group")
    elif close_sparse:
        flags.append("close_group_sparse")
    if long_branch:
        flags.append("query_long_branch")
    if query_branch_ratio != math.inf and query_branch_ratio >= args.long_branch_ratio:
        flags.append("query_branch_high_ratio")
    if orf_bad:
        flags.append("query_orf_internal_stops_or_frame_uncertain")
    elif orf_minor:
        flags.append("query_orf_minor_stop")
    if local_mixed:
        flags.append("donor_recipient_mixed_local_split")
    if nearest_mixed:
        flags.append("donor_recipient_mixed_nearest10")
    if sup_level == "low":
        flags.append("local_split_low_support")
    elif sup_level == "none":
        flags.append("no_small_query_split_support")
    if donor_count < args.min_donor_context:
        flags.append("few_donor_sequences")
    if counts.get("unknown", 0) + counts.get("other", 0) >= args.many_unknown_cutoff:
        flags.append("many_unknown_or_other_sequences")

    # Primary classification. Conservative and multi-tiered.
    # v2 logic: DO NOT call vertical inheritance when CLOSE sampling is absent/sparse.
    # Biological rationale: without close/outside-recipient eukaryotic homologs,
    # QUERY clustering with SELF can represent vertical inheritance, but it can also
    # represent an older transfer into the recipient-lineage ancestor followed by
    # retention. Therefore, no-close trees are retained as possible/ambiguous HGT
    # unless long-branch/ORF artifacts dominate.
    primary_class = "unclassified"
    evidence_tier = "unresolved"
    hgt_direction = "unresolved"
    reason = []

    no_close = close_count == 0
    very_sparse_close = close_count <= args.few_close

    if donor_count == 0:
        primary_class = "no_donor_context"
        evidence_tier = "not_evaluable"
        hgt_direction = "not_evaluable"
        reason.append("No donor/out-group homologs are present in the tree.")

    elif orf_bad and (long_branch or donor_nearer or nearest_mixed):
        primary_class = "translation_frame_or_ORF_artifact_suspected"
        evidence_tier = "not_reliable"
        hgt_direction = "artifact_possible"
        reason.append("Query ORF has internal stops/frame uncertainty and tree placement may be affected.")

    elif no_close:
        # No CLOSE group: avoid vertical-inheritance calls.
        # This is the key change from v1.
        if long_branch:
            primary_class = "long_branch_no_close_abnormal"
            evidence_tier = "not_reliable_or_manual_review"
            hgt_direction = "artifact_possible"
            reason.append("No CLOSE homologs are available and QUERY is an abnormal long branch; HGT cannot be distinguished from long-branch/artifact.")
        elif orf_bad:
            primary_class = "no_close_ORF_uncertain_possible_HGT"
            evidence_tier = "weak_possible_HGT_but_ORF_uncertain"
            hgt_direction = "HGT_possible_but_ORF_uncertain"
            reason.append("No CLOSE homologs are available; QUERY is not long-branched but ORF/frame quality is uncertain.")
        elif donor_only_context:
            primary_class = "donor_only_context_possible_HGT"
            evidence_tier = "weak_possible_HGT_close_absent"
            hgt_direction = "HGT_possible_but_no_close_contrast"
            reason.append("Only donor/out-group context is available; in the absence of CLOSE homologs this is retained as possible HGT if no artifact is detected.")
        elif donor_nearer or near10_donor >= max(2, near10_recipient + 1):
            primary_class = "possible_HGT_no_close_donor_leaning"
            evidence_tier = "weak_possible_HGT_close_absent"
            hgt_direction = "HGT_leaning_but_no_close_contrast"
            reason.append("No CLOSE homologs are available and QUERY is donor-near; retained as possible HGT but not strong due to missing close contrast.")
        elif recipient_nearer or recipient_only_context or near10_recipient > 0:
            primary_class = "recipient_lineage_like_no_close_possible_ancestral_HGT"
            evidence_tier = "ambiguous_possible_ancestral_HGT_or_vertical"
            hgt_direction = "possible_ancestral_HGT_or_vertical_unresolved"
            reason.append("No CLOSE homologs are available; QUERY is close to SELF/recipient-lineage homologs, which is compatible with vertical inheritance but also with ancestral HGT followed by retention.")
        elif nearest_mixed or local_mixed:
            primary_class = "no_close_donor_recipient_mixed_possible_HGT"
            evidence_tier = "ambiguous_possible_HGT_close_absent"
            hgt_direction = "unresolved_possible_HGT"
            reason.append("No CLOSE homologs are available and donor/recipient leaves are locally mixed; retained for manual review as possible HGT.")
        else:
            primary_class = "possible_HGT_no_close_unresolved"
            evidence_tier = "weak_possible_HGT_close_absent"
            hgt_direction = "HGT_possible_but_unresolved_no_close"
            reason.append("No CLOSE homologs are available; no long-branch/ORF artifact is detected, but topology lacks close contrast.")

    elif very_sparse_close:
        # Sparse CLOSE group: also avoid strong vertical calls.
        if long_branch:
            primary_class = "long_branch_sparse_close_abnormal"
            evidence_tier = "not_reliable_or_manual_review"
            hgt_direction = "artifact_possible"
            reason.append("QUERY is a long branch and close-group sampling is sparse.")
        elif donor_nearer and near10_donor >= args.min_near_donor_for_hgt:
            primary_class = "possible_HGT_close_sparse"
            evidence_tier = "weak_possible_HGT_close_sparse"
            hgt_direction = "HGT_leaning_but_close_sparse"
            reason.append("QUERY is donor-near, but close-group sampling is sparse; HGT remains possible but not strong.")
        elif recipient_nearer or recipient_only_context:
            primary_class = "recipient_lineage_like_close_sparse_possible_ancestral_HGT"
            evidence_tier = "ambiguous_possible_ancestral_HGT_or_vertical_close_sparse"
            hgt_direction = "possible_ancestral_HGT_or_vertical_close_sparse"
            reason.append("QUERY is closest to available SELF/CLOSE homologs, but close sampling is sparse, so this is not called strong vertical inheritance.")
        elif nearest_mixed or local_mixed:
            primary_class = "sparse_close_mixed_unresolved"
            evidence_tier = "ambiguous_close_sparse"
            hgt_direction = "unresolved_conflict_close_sparse"
            reason.append("Sparse close-group context and mixed donor/recipient topology.")
        else:
            primary_class = "sparse_close_unresolved"
            evidence_tier = "unresolved_close_sparse"
            hgt_direction = "unresolved_close_sparse"
            reason.append("Close-group sampling is sparse and no reliable donor/recipient signal is available.")

    else:
        # close-rich / informative context: can separate HGT vs vertical-like.
        if long_branch:
            primary_class = "long_branch_artifact_possible"
            evidence_tier = "manual_review"
            hgt_direction = "artifact_possible"
            reason.append("QUERY terminal branch is much longer than typical terminal branches.")
        elif donor_nearer and near10_donor >= max(args.min_near_donor_for_hgt, near10_recipient + 1) and near5_recipient <= args.max_near5_recipient_for_hgt:
            if sup_level == "strong" and not orf_minor:
                primary_class = "reliable_HGT"
                evidence_tier = "strong_HGT"
            else:
                primary_class = "probable_HGT"
                evidence_tier = "probable_HGT_support_limited"
            hgt_direction = "HGT_leaning"
            reason.append("With informative CLOSE context, QUERY is closer to donor/out-group sequences than to recipient SELF/CLOSE homologs.")
        elif recipient_nearer and near5_recipient >= args.min_near5_recipient_for_vertical and near5_donor == 0:
            if sup_level in {"strong", "moderate", "none"} and not local_mixed:
                primary_class = "vertical_inheritance_likely"
                evidence_tier = "vertical_like"
            else:
                primary_class = "vertical_inheritance_likely_low_support"
                evidence_tier = "vertical_like_low_support"
            hgt_direction = "vertical_leaning"
            reason.append("With informative CLOSE/recipient context, QUERY is closest to SELF/CLOSE recipient homologs and not locally donor-dominated.")
        elif local_mixed or nearest_mixed:
            if donor_nearer and near10_donor > near10_recipient:
                primary_class = "conflicting_topology_HGT_leaning"
                evidence_tier = "ambiguous_HGT_leaning"
                hgt_direction = "HGT_leaning_but_conflicted"
            elif recipient_nearer and near10_recipient >= near10_donor:
                primary_class = "conflicting_topology_vertical_leaning"
                evidence_tier = "ambiguous_vertical_leaning"
                hgt_direction = "vertical_leaning_but_conflicted"
            else:
                primary_class = "donor_recipient_mixed_unresolved"
                evidence_tier = "ambiguous_conflicting"
                hgt_direction = "unresolved_conflict"
            reason.append("Donor and recipient sequences are mixed in the QUERY neighborhood or small split.")
        elif donor_nearer:
            primary_class = "possible_HGT_low_local_support"
            evidence_tier = "weak_possible_HGT"
            hgt_direction = "HGT_leaning"
            reason.append("Donor sequences are nearer, but local dominance/support criteria are not strong.")
        elif recipient_nearer or recipient_only_context:
            primary_class = "vertical_inheritance_possible"
            evidence_tier = "vertical_like_weak"
            hgt_direction = "vertical_leaning"
            reason.append("Recipient SELF/CLOSE sequences are nearer, but the topology is not clean enough for a strong call.")
        else:
            primary_class = "close_rich_unresolved"
            evidence_tier = "unresolved_with_close_context"
            hgt_direction = "unresolved"
            reason.append("Close/recipient homologs are present, but donor/recipient distances and local topology do not give a clear direction.")

    if orf_minor and primary_class not in {"translation_frame_or_ORF_artifact_suspected"}:
        flags.append("evidence_downgraded_for_minor_ORF_issue")
    if orf_bad and primary_class not in {"translation_frame_or_ORF_artifact_suspected"}:
        flags.append("do_not_use_as_strong_HGT_due_to_ORF")

    return {
        "hgt_id": hgt_id,
        "primary_class": primary_class,
        "evidence_tier": evidence_tier,
        "hgt_direction": hgt_direction,
        "reason": " ".join(reason),
        "flags": ";".join(sorted(set(flags))),
        "query_name": query_name,
        "leaf_n": len(lf),
        "self_n": self_count,
        "close_n": close_count,
        "recipient_n": recipient_count,
        "donor_n": donor_count,
        "unknown_other_n": counts.get("unknown", 0) + counts.get("other", 0),
        "nearest1_cat": nearest1_cat,
        "nearest1_dist": fmt_num(nearest1_dist),
        "nearest1_name": nearest1_name,
        "nearest5_recipient_n": near5_recipient,
        "nearest5_donor_n": near5_donor,
        "nearest10_recipient_n": near10_recipient,
        "nearest10_donor_n": near10_donor,
        "min_dist_self": fmt_num(d_self),
        "min_dist_close": fmt_num(d_close),
        "min_dist_recipient": fmt_num(d_recipient),
        "min_dist_donor": fmt_num(d_donor),
        "min_dist_unknown": fmt_num(d_unknown),
        "query_branch_length": fmt_num(query_branch),
        "median_terminal_branch_length": fmt_num(median_terminal),
        "query_branch_ratio_to_median": fmt_num(query_branch_ratio),
        "long_branch": str(long_branch),
        "support_level_query_split": sup_level,
        "query_split_size": split.get("query_split_size"),
        "query_split_self": split.get("query_split_self"),
        "query_split_close": split.get("query_split_close"),
        "query_split_donor": split.get("query_split_donor"),
        "query_split_unknown": split.get("query_split_unknown"),
        "query_split_sh_alrt": fmt_num(split.get("query_split_sh_alrt")),
        "query_split_ufboot": fmt_num(split.get("query_split_ufboot")),
        "orf_known": str(orf_known),
        "orf_status": orf_status,
        "frame_status": frame_status,
        "best_strand": best_strand,
        "best_frame_0based": best_frame,
        "stop_codons_to_X": stop_codons if stop_codons >= 0 else "NA",
        "ambiguous_codons_to_X": amb_codons if amb_codons >= 0 else "NA",
        "treefile": str(treefile),
    }


def fmt_num(x: object) -> str:
    try:
        if x is None:
            return "NA"
        xf = float(x)
        if not math.isfinite(xf):
            return "NA"
        return f"{xf:.8g}"
    except Exception:
        if x in {"NA", ""}:
            return "NA"
        return str(x)


# -----------------------------
# Missing tree diagnostics
# -----------------------------

def write_no_tree_report(jobs_dir: Path, tree_ids: Set[str], out_path: Path) -> None:
    rows = []
    for d in sorted(p for p in jobs_dir.iterdir() if p.is_dir()):
        hgt_id = d.name
        if hgt_id in tree_ids:
            continue
        fa = d / "input.aa.fa"
        aln = d / "aln.fa"
        trim = d / "aln.trim.fa"
        status = d / "tree.status"
        if not fa.exists() or fa.stat().st_size == 0:
            reason = "missing_input"
        elif not aln.exists() or aln.stat().st_size == 0:
            reason = status.read_text().strip() if status.exists() else "not_started_or_no_status"
        elif not trim.exists() or trim.stat().st_size == 0:
            reason = status.read_text().strip() if status.exists() else "mafft_done_trim_not_done"
        else:
            reason = status.read_text().strip() if status.exists() else "trim_done_iqtree_not_done"
        rows.append({"hgt_id": hgt_id, "reason": reason})
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["hgt_id", "reason"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Classify HGT candidate phylogenies with no-close-aware conservative logic; no-CLOSE trees are not called vertical but retained as possible/ambiguous HGT unless artifact-prone."
    )
    ap.add_argument("--outdir", default="hgt_phylo_out", help="Pipeline output directory containing jobs/. Default: hgt_phylo_out")
    ap.add_argument("--jobs-dir", default=None, help="Override jobs directory. Default: <outdir>/jobs")
    ap.add_argument("--orf-summary", default=None, help="Optional nucleotide ORF/frame summary TSV. Default: <outdir>/nucleotide_orf_priority_frame_summary.tsv if present")
    ap.add_argument("--output-dir", default=None, help="Output directory. Default: <outdir>/tree_classification")

    # Sampling context thresholds
    ap.add_argument("--min-close-rich", type=int, default=5, help="close_n >= this is treated as close-rich context. Default: 5")
    ap.add_argument("--few-close", type=int, default=2, help="close_n <= this is treated as sparse close context. Default: 2")
    ap.add_argument("--min-recipient-rich", type=int, default=5, help="self_n+close_n >= this gives recipient-rich context. Default: 5")
    ap.add_argument("--min-donor-context", type=int, default=2, help="Flag trees with fewer donor sequences than this. Default: 2")
    ap.add_argument("--many-unknown-cutoff", type=int, default=20, help="Flag trees with >= this unknown/other leaves. Default: 20")

    # Neighborhood / distance thresholds
    ap.add_argument("--nearest-k", type=int, default=10, help="Number of nearest leaves summarized around QUERY. Default: 10")
    ap.add_argument("--distance-ratio", type=float, default=0.80, help="A side is nearer if min distance < other side * ratio. Default: 0.80")
    ap.add_argument("--min-near-donor-for-hgt", type=int, default=3, help="Minimum donor leaves among nearest10 for HGT leaning calls. Default: 3")
    ap.add_argument("--max-near5-recipient-for-hgt", type=int, default=1, help="Max recipient leaves among nearest5 for reliable/probable HGT. Default: 1")
    ap.add_argument("--min-near5-recipient-for-vertical", type=int, default=3, help="Minimum recipient leaves among nearest5 for vertical-like calls. Default: 3")
    ap.add_argument("--max-local-split", type=int, default=25, help="Max query-side split size used to inspect local topology. Default: 25")

    # Long branch thresholds
    ap.add_argument("--long-branch-ratio", type=float, default=10.0, help="QUERY terminal branch / median terminal branch ratio threshold. Default: 10")
    ap.add_argument("--long-branch-abs", type=float, default=2.0, help="Absolute query branch threshold combined with ratio. Default: 2.0")
    ap.add_argument("--very-long-branch-abs", type=float, default=5.0, help="Always flag query branch >= this. Default: 5.0")

    # Support thresholds
    ap.add_argument("--sh-strong", type=float, default=80.0, help="Strong SH-aLRT threshold. Default: 80")
    ap.add_argument("--uf-strong", type=float, default=95.0, help="Strong UFBoot threshold. Default: 95")
    ap.add_argument("--sh-moderate", type=float, default=70.0, help="Moderate SH-aLRT threshold. Default: 70")
    ap.add_argument("--uf-moderate", type=float, default=80.0, help="Moderate UFBoot threshold. Default: 80")

    # ORF thresholds
    ap.add_argument("--minor-stop-cutoff", type=int, default=1, help="stop_codons_to_X >= this is a minor ORF issue. Default: 1")
    ap.add_argument("--bad-stop-cutoff", type=int, default=2, help="stop_codons_to_X >= this is a serious ORF/frame issue. Default: 2")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    jobs_dir = Path(args.jobs_dir) if args.jobs_dir else outdir / "jobs"
    output_dir = Path(args.output_dir) if args.output_dir else outdir / "tree_classification"
    output_dir.mkdir(parents=True, exist_ok=True)

    orf_path = Path(args.orf_summary) if args.orf_summary else outdir / "nucleotide_orf_priority_frame_summary.tsv"
    orf_by_id = read_tsv_by_key(orf_path if orf_path.exists() else None, "hgt_id")

    treefiles = sorted(jobs_dir.glob("*/iqtree.treefile"))
    rows: List[Dict[str, object]] = []
    errors: List[Dict[str, str]] = []

    print(f"[INFO] jobs_dir: {jobs_dir}")
    print(f"[INFO] treefiles: {len(treefiles)}")
    print(f"[INFO] ORF summary: {orf_path if orf_path.exists() else 'not_found'} rows={len(orf_by_id)}")

    for i, tf in enumerate(treefiles, start=1):
        hgt_id = parse_hgt_id_from_treefile(tf)
        try:
            row = summarize_tree(tf, args, orf_by_id.get(hgt_id, {}))
            rows.append(row)
        except Exception as e:
            errors.append({"hgt_id": hgt_id, "treefile": str(tf), "error": repr(e)})
        if i % 1000 == 0:
            print(f"[INFO] classified {i}/{len(treefiles)} trees; errors={len(errors)}", flush=True)

    out_tsv = output_dir / "tree_classification.tsv"
    if rows:
        fieldnames = list(rows[0].keys())
        with out_tsv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    print(f"[INFO] wrote: {out_tsv} rows={len(rows)}")

    err_tsv = output_dir / "tree_classification_errors.tsv"
    with err_tsv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["hgt_id", "treefile", "error"], delimiter="\t")
        w.writeheader()
        w.writerows(errors)
    print(f"[INFO] wrote: {err_tsv} errors={len(errors)}")

    # Summary tables
    def write_counter(col: str, outname: str):
        c = Counter(str(r.get(col, "")) for r in rows)
        with (output_dir / outname).open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow([col, "n"])
            for k, v in c.most_common():
                w.writerow([k, v])
        return c

    class_counts = write_counter("primary_class", "summary_by_primary_class.tsv")
    write_counter("evidence_tier", "summary_by_evidence_tier.tsv")
    write_counter("hgt_direction", "summary_by_hgt_direction.tsv")

    # Per-class ID lists
    class_dir = output_dir / "class_ids"
    class_dir.mkdir(exist_ok=True)
    by_class: Dict[str, List[str]] = defaultdict(list)
    for r in rows:
        by_class[str(r["primary_class"])].append(str(r["hgt_id"]))
    for cls, ids in by_class.items():
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", cls)
        (class_dir / f"{safe}.ids.txt").write_text("\n".join(ids) + "\n")

    # Manual review subset
    manual_keywords = ["ambiguous", "conflicting", "unresolved", "artifact", "long_branch", "sparse", "ORF", "mixed"]
    manual_rows = []
    for r in rows:
        text = (str(r.get("primary_class", "")) + ";" + str(r.get("evidence_tier", "")) + ";" + str(r.get("flags", ""))).lower()
        if any(k.lower() in text for k in manual_keywords):
            manual_rows.append(r)
    manual_tsv = output_dir / "manual_review.tsv"
    if rows:
        with manual_tsv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(manual_rows)
    print(f"[INFO] wrote: {manual_tsv} rows={len(manual_rows)}")

    # Missing tree diagnostics
    no_tree_tsv = output_dir / "no_tree_report.tsv"
    tree_ids = {str(r["hgt_id"]) for r in rows}
    if jobs_dir.exists():
        write_no_tree_report(jobs_dir, tree_ids, no_tree_tsv)
        print(f"[INFO] wrote: {no_tree_tsv}")

    # Console summary
    print("\n[SUMMARY primary_class]")
    for k, v in class_counts.most_common():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
