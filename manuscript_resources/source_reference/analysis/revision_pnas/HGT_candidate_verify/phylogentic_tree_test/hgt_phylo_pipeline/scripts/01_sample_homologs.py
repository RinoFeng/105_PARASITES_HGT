#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import re
from typing import Dict, List, Tuple
import pandas as pd

from hgt_tree_utils import (
    read_taxonomy_file,
    read_genome_taxid_map,
    get_self_taxid_for_target,
    classify_hit,
    taxonomy_rank,
    species_name,
    genus_name,
)


def align_paths(align_dir: str, genome: str):
    return (
        os.path.join(align_dir, f'{genome}.out'),
        os.path.join(align_dir, f'{genome}.acc2taxid.out')
    )




_COORD_RE = re.compile(r'^(?P<seqid>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)(?:\([+-]\))?$')

def normalize_query_id(q: str) -> str:
    """Normalize nucleotide region IDs for matching.

    Examples that should map together:
      CM000442.1:40422-40783
      CM000442.1:40783-40422
      CM000442.1:40422-40783(+)
    Protein accessions are returned unchanged except whitespace stripping.
    """
    q = str(q).strip()
    q = q.split()[0]
    q = q.replace('..', '-')
    q = re.sub(r'\([+-]\)$', '', q)
    m = _COORD_RE.match(q)
    if not m:
        return q
    seqid = m.group('seqid')
    a = int(m.group('start'))
    b = int(m.group('end'))
    lo, hi = sorted((a, b))
    return f'{seqid}:{lo}-{hi}'

def parse_region_id(q: str):
    qn = normalize_query_id(q)
    m = _COORD_RE.match(qn)
    if not m:
        return None
    seqid = m.group('seqid')
    start = int(m.group('start'))
    end = int(m.group('end'))
    lo, hi = sorted((start, end))
    return seqid, lo, hi

def interval_overlap_fraction(a1: int, a2: int, b1: int, b2: int) -> float:
    a1, a2 = sorted((int(a1), int(a2)))
    b1, b2 = sorted((int(b1), int(b2)))
    ov = max(0, min(a2, b2) - max(a1, b1) + 1)
    denom = max(1, min(a2-a1+1, b2-b1+1))
    return ov / denom

def load_hit_table(out_file: str, acc2taxid_file: str) -> pd.DataFrame:
    """Load DIAMOND tabular hits and taxids, then merge.

    .out columns:
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    .acc2taxid.out columns:
      qseqid qtaxid sseqid staxid pident length bitscore

    Notes:
      - protein .out in the user's files may store pident as fraction, e.g. 0.964.
      - nucleotide/blastx .out may store pident as percent, e.g. 99.2.
      This script normalizes pident to percent units before filtering.
    """
    if not os.path.exists(out_file):
        raise FileNotFoundError(out_file)
    if not os.path.exists(acc2taxid_file):
        raise FileNotFoundError(acc2taxid_file)

    out_cols = ['qseqid','sseqid','pident','aln_len','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    tax_cols = ['qseqid','hit_qtaxid','sseqid','target_taxid','pident_tax','aln_len_tax','bitscore_tax']

    out = pd.read_csv(out_file, sep='\t', names=out_cols, dtype=str)
    tax = pd.read_csv(acc2taxid_file, sep='\t', names=tax_cols, dtype=str)
    tax = tax.drop_duplicates(['qseqid','sseqid'], keep='first')

    df = out.merge(tax[['qseqid','sseqid','hit_qtaxid','target_taxid']], on=['qseqid','sseqid'], how='left')
    for c in ['pident','aln_len','qstart','qend','sstart','send','bitscore']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['evalue_num'] = pd.to_numeric(df['evalue'], errors='coerce')

    # Normalize pident to percent. Protein DIAMOND output in the user's example is 0.964.
    max_pid = df['pident'].max(skipna=True)
    if pd.notna(max_pid) and max_pid <= 1.5:
        df['pident'] = df['pident'] * 100.0

    # Normalized query IDs are critical for nucleotide/blastx rows, where region IDs
    # may differ by strand suffix or coordinate order. This does not change protein IDs.
    df['qseqid_norm'] = df['qseqid'].astype(str).map(normalize_query_id)

    # Parsed nucleotide query region fields; protein rows remain NA.
    parsed = df['qseqid_norm'].map(parse_region_id)
    df['qseqid_contig'] = parsed.map(lambda x: x[0] if x else '')
    df['qseqid_start'] = pd.to_numeric(parsed.map(lambda x: x[1] if x else None), errors='coerce')
    df['qseqid_end'] = pd.to_numeric(parsed.map(lambda x: x[2] if x else None), errors='coerce')

    return df


def choose_diverse(df: pd.DataFrame, max_n: int) -> pd.DataFrame:
    """Choose taxonomically diverse high-scoring hits.

    The goal is a tree-informative set, not a top-hit pile:
      1. one representative per genus/species/taxid if possible
      2. fill remaining by best score while avoiding exact duplicate accession
    """
    if df.empty or max_n <= 0:
        return df.head(0)
    df = df.sort_values(['bitscore','aln_len','pident'], ascending=[False, False, False]).copy()
    df = df.drop_duplicates(['sseqid'], keep='first')

    chosen_idx: List[int] = []
    used_taxid, used_species, used_genus = set(), set(), set()

    for idx, row in df.iterrows():
        if len(chosen_idx) >= max_n:
            break
        tx = str(row.get('target_taxid',''))
        sp = str(row.get('species',''))
        gn = str(row.get('genus',''))
        if tx not in used_taxid and sp not in used_species and gn not in used_genus:
            chosen_idx.append(idx)
            used_taxid.add(tx); used_species.add(sp); used_genus.add(gn)

    if len(chosen_idx) < max_n:
        for idx, row in df.iterrows():
            if len(chosen_idx) >= max_n:
                break
            if idx in chosen_idx:
                continue
            tx = str(row.get('target_taxid',''))
            sp = str(row.get('species',''))
            if tx not in used_taxid and sp not in used_species:
                chosen_idx.append(idx)
                used_taxid.add(tx); used_species.add(sp)

    if len(chosen_idx) < max_n:
        for idx, _ in df.iterrows():
            if len(chosen_idx) >= max_n:
                break
            if idx not in chosen_idx:
                chosen_idx.append(idx)

    return df.loc[chosen_idx].copy()




def get_query_hits(hits: pd.DataFrame,
                   hgt_row: pd.Series,
                   min_overlap: float = 0.50,
                   coord_tolerance: int = 2) -> Tuple[pd.DataFrame, str, str]:
    """Return hits for one HGT target using robust query-ID matching.

    Matching order:
      1. exact Accession match
      2. normalized Accession match
      3. exact/normalized `query` column match
      4. nucleotide coordinate-tolerance match on the same contig
         Example: curated CM000444.2:175849-176494 matches blast qseqid
         CM000444.2:175848-176494 when coord_tolerance >= 1.
      5. coordinate-overlap fallback on the same contig

    The tolerance step is intentionally more specific than generic overlap and is
    useful for 0/1-based coordinate differences introduced during region export.
    """
    q = str(hgt_row.get('Accession', '')).strip()
    if not q or q.lower() == 'nan':
        return hits.head(0).copy(), '', 'missing_accession'

    qhits = hits[hits['qseqid'].astype(str) == q].copy()
    if not qhits.empty:
        return qhits, q, 'exact_accession'

    q_norm = normalize_query_id(q)
    if 'qseqid_norm' in hits.columns:
        qhits = hits[hits['qseqid_norm'] == q_norm].copy()
        if not qhits.empty:
            return qhits, qhits['qseqid'].iloc[0], 'normalized_accession'

    q2 = str(hgt_row.get('query', '')).strip()
    if q2 and q2.lower() != 'nan' and q2 != q:
        qhits = hits[hits['qseqid'].astype(str) == q2].copy()
        if not qhits.empty:
            return qhits, q2, 'exact_query_column'
        q2_norm = normalize_query_id(q2)
        if 'qseqid_norm' in hits.columns:
            qhits = hits[hits['qseqid_norm'] == q2_norm].copy()
            if not qhits.empty:
                return qhits, qhits['qseqid'].iloc[0], 'normalized_query_column'

    # Coordinate-based fallback for Nucleotide region IDs.
    # First try a strict endpoint tolerance, then a broader overlap fallback.
    region = parse_region_id(q)
    if region and {'qseqid_contig','qseqid_start','qseqid_end'}.issubset(hits.columns):
        contig, start, end = region
        sub0 = hits[hits['qseqid_contig'] == contig].copy()

        if not sub0.empty:
            tol = int(coord_tolerance or 0)
            sub = sub0[
                (sub0['qseqid_start'].notna()) &
                (sub0['qseqid_end'].notna()) &
                ((sub0['qseqid_start'] - start).abs() <= tol) &
                ((sub0['qseqid_end'] - end).abs() <= tol)
            ].copy()
            if not sub.empty:
                # If several qseqids satisfy tolerance, pick the one with the
                # smallest endpoint error, then the strongest total bitscore.
                sub['_coord_delta'] = (sub['qseqid_start'] - start).abs() + (sub['qseqid_end'] - end).abs()
                q_rank = (
                    sub.groupby('qseqid')
                       .agg(coord_delta=('_coord_delta','min'), bitscore_sum=('bitscore','sum'))
                       .sort_values(['coord_delta','bitscore_sum'], ascending=[True, False])
                )
                best_q = q_rank.index[0]
                qhits = sub[sub['qseqid'] == best_q].drop(columns=['_coord_delta'], errors='ignore').copy()
                return qhits, best_q, f'coordinate_tolerance_{tol}bp_accession'

            ov = sub0.apply(
                lambda r: interval_overlap_fraction(start, end, r['qseqid_start'], r['qseqid_end'])
                if pd.notna(r['qseqid_start']) and pd.notna(r['qseqid_end']) else 0.0,
                axis=1
            )
            sub = sub0[ov >= min_overlap].copy()
            if not sub.empty:
                sub['_query_overlap'] = ov.loc[sub.index]
                q_rank = (
                    sub.groupby('qseqid')
                       .agg(query_overlap=('_query_overlap','max'), bitscore_sum=('bitscore','sum'))
                       .sort_values(['query_overlap','bitscore_sum'], ascending=[False, False])
                )
                best_q = q_rank.index[0]
                qhits = sub[sub['qseqid'] == best_q].drop(columns=['_query_overlap'], errors='ignore').copy()
                return qhits, best_q, 'coordinate_overlap_accession'

    return hits.head(0).copy(), q, 'no_query_hit'

def to_int(x, default=0) -> int:
    try:
        if pd.isna(x):
            return default
        return int(float(str(x)))
    except Exception:
        return default


def sample_for_target(hgt_row: pd.Series,
                      hits: pd.DataFrame,
                      tax_info: Dict[str, Dict[str,str]],
                      genome_taxid_map: Dict[str, Dict[str, Dict[str,str]]],
                      args) -> Tuple[pd.DataFrame, Dict[str, object]]:
    q = str(hgt_row['Accession'])
    hgt_id = str(hgt_row['hgt_id'])
    qhits, matched_query_id, query_match_mode = get_query_hits(hits, hgt_row, min_overlap=args.query_overlap, coord_tolerance=args.query_coord_tolerance)

    meta = {
        'hgt_id': hgt_id,
        'Data_type': str(hgt_row.get('Data_type','')),
        'Genome': str(hgt_row.get('Genome','')),
        'Accession': q,
        'matched_query_id': matched_query_id,
        'query_match_mode': query_match_mode,
        'tree_priority': str(hgt_row.get('tree_priority','')),
        'raw_hit_n': int(len(qhits)),
        'passed_filter_n': 0,
        'self_taxid': '',
        'self_taxid_source': '',
        'self_level': '',
        'self_n': 0,
        'close_n': 0,
        'out_host_n': 0,
        'out_other_n': 0,
        'unknown_n': 0,
        'selected_n': 0,
        'status': 'init',
    }

    if qhits.empty:
        meta['status'] = query_match_mode if query_match_mode else 'no_query_hit'
        return qhits, meta

    qhits = qhits[(qhits['aln_len'].fillna(0) >= args.min_align_len) & (qhits['pident'].fillna(0) >= args.min_pident)].copy()
    meta['passed_filter_n'] = int(len(qhits))
    if qhits.empty:
        meta['status'] = 'no_hit_after_len_identity_filter'
        return qhits, meta

    # CRITICAL FIX: recipient/self taxid comes from 105_name2accfil_withTaxid.txt,
    # matching HGT_verify_parasite.py usage:
    #   python HGT_verify_parasite.py taxonomy.txt genome.acc2taxid.out "$taxid"
    self_taxid, self_taxid_source = get_self_taxid_for_target(hgt_row, genome_taxid_map)
    meta['self_taxid'] = self_taxid
    meta['self_taxid_source'] = self_taxid_source

    if not self_taxid or str(self_taxid) not in tax_info:
        meta['status'] = 'missing_self_taxid_or_self_taxonomy'
        # Do not silently sample; this would make trees biologically meaningless.
        return qhits.head(0), meta

    self_info = tax_info.get(str(self_taxid), {'detail':'','short':''})
    try:
        if ';Metazoa;' in self_info.get('detail',''):
            meta['self_level'] = self_info.get('short','').split(';')[1]
        else:
            meta['self_level'] = self_info.get('detail','').split(';')[2]
    except Exception:
        meta['self_level'] = ''

    groups, details, shorts, species, genus, rank2 = [], [], [], [], [], []
    for tx in qhits['target_taxid'].astype(str):
        g, detail, short = classify_hit(tax_info, self_taxid, tx)
        groups.append(g)
        details.append(detail)
        shorts.append(short)
        species.append(species_name(short))
        genus.append(genus_name(short))
        rank2.append(taxonomy_rank(short, 2))

    qhits['hit_group'] = groups
    qhits['taxonomy_detail'] = details
    qhits['taxonomy_short'] = shorts
    qhits['species'] = species
    qhits['genus'] = genus
    qhits['rank2'] = rank2
    qhits['hgt_id'] = hgt_id
    qhits['data_type'] = str(hgt_row['Data_type'])
    qhits['tree_priority'] = str(hgt_row.get('tree_priority',''))
    qhits['phylo_testability'] = str(hgt_row.get('phylo_testability',''))
    qhits['main_analysis_set'] = str(hgt_row.get('main_analysis_set',''))
    qhits['recipient_self_taxid'] = self_taxid
    qhits['recipient_self_level'] = meta['self_level']

    counts = qhits['hit_group'].value_counts(dropna=False).to_dict()
    meta['self_n'] = int(counts.get('self', 0))
    meta['close_n'] = int(counts.get('close', 0))
    meta['out_host_n'] = int(counts.get('out_host', 0))
    meta['out_other_n'] = int(counts.get('out_other', 0))
    meta['unknown_n'] = int(counts.get('unknown', 0))

    self_df = qhits[qhits.hit_group == 'self']
    close_df = qhits[qhits.hit_group == 'close']
    out_df = qhits[qhits.hit_group.isin(['out_other','out_host'])]

    close_n_from_curated = to_int(hgt_row.get('Close_num', 0), 0)
    tree_priority = str(hgt_row.get('tree_priority',''))

    max_self = args.max_self
    max_close = args.max_close
    max_out = args.max_out

    # Dynamic sampling. For close-absent candidates, only sample self + donor/out groups.
    if 'P1' in tree_priority:
        max_close = max(args.max_close, 35)
        max_out = min(max(args.max_out, 70), max(0, args.max_total - max_self - max_close))
    elif close_n_from_curated == 0:
        max_close = 0
        max_self = max(args.max_self, 15)
        max_out = min(max(args.max_out, 90), max(0, args.max_total - max_self))
    elif close_n_from_curated <= 2:
        max_close = min(args.max_close, 12)
        max_out = min(args.max_out, 70)

    chosen = []
    if not self_df.empty:
        chosen.append(choose_diverse(self_df, max_self))
    if max_close > 0 and not close_df.empty:
        chosen.append(choose_diverse(close_df, max_close))
    if not out_df.empty:
        chosen.append(choose_diverse(out_df, max_out))

    selected = pd.concat(chosen, ignore_index=True) if chosen else qhits.head(0)

    # Optional diagnostic fallback, disabled by default. Useful only to inspect query hits.
    if selected.empty and args.fallback_top_hits > 0:
        selected = qhits.sort_values(['bitscore','aln_len','pident'], ascending=[False,False,False]).head(args.fallback_top_hits).copy()
        selected['hit_group'] = selected['hit_group'].replace({'unknown':'fallback_unknown'})
        meta['status'] = 'fallback_top_hits_used'

    if not selected.empty:
        group_order = {'self':0, 'close':1, 'out_other':2, 'out_host':3, 'unknown':4, 'fallback_unknown':5}
        selected['_group_order'] = selected['hit_group'].map(group_order).fillna(9)
        selected = selected.sort_values(['_group_order','bitscore','aln_len'], ascending=[True, False, False])
        if len(selected) > args.max_total:
            s1 = selected[selected.hit_group == 'self'].head(max_self)
            s2 = selected[selected.hit_group == 'close'].head(max_close)
            remaining = args.max_total - len(s1) - len(s2)
            s3 = selected[selected.hit_group.isin(['out_other','out_host'])].head(max(0, remaining))
            selected = pd.concat([s1, s2, s3], ignore_index=True)
        selected = selected.drop(columns=['_group_order'], errors='ignore')
        selected['tree_role'] = selected['hit_group'].replace({
            'self':'SELF',
            'close':'CLOSE',
            'out_other':'DONOR_OUT',
            'out_host':'HOST_OUT',
            'fallback_unknown':'FALLBACK_UNKNOWN'
        })

    meta['selected_n'] = int(len(selected))
    if meta['status'] == 'init':
        meta['status'] = 'ok' if len(selected) > 0 else 'no_selected_hit'
    return selected, meta


def main():
    ap = argparse.ArgumentParser(description='Sample self/close/out homologs for HGT gene trees using HGT_verify_parasite.py taxonomy rules.')
    ap.add_argument('--targets', required=True)
    ap.add_argument('--genome-taxid-map', required=True, help='105_name2accfil_withTaxid.txt: species, genome FASTA, taxid')
    ap.add_argument('--protein-align-dir', required=True)
    ap.add_argument('--protein-taxonomy', required=True)
    ap.add_argument('--nt-align-dir', required=True)
    ap.add_argument('--nt-taxonomy', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--max-self', type=int, default=12)
    ap.add_argument('--max-close', type=int, default=30)
    ap.add_argument('--max-out', type=int, default=70)
    ap.add_argument('--max-remote-euk', type=int, default=10)  # kept for CLI compatibility; not used by original-rule sampler
    ap.add_argument('--max-total', type=int, default=130)
    ap.add_argument('--min-align-len', type=int, default=80)
    ap.add_argument('--min-pident', type=float, default=20)
    ap.add_argument('--fallback-top-hits', type=int, default=0, help='Diagnostic only. If >0, sample top N hits when grouping returns nothing.')
    ap.add_argument('--query-overlap', type=float, default=0.50, help='For nucleotide region IDs, allow qseqid matching by same contig and interval overlap fraction.')
    ap.add_argument('--query-coord-tolerance', type=int, default=2, help='For nucleotide region IDs, allow same-contig start/end differences up to this many bp before overlap fallback.')
    ap.add_argument('--only-data-type', choices=['Protein','Nucleotide'], default='', help='Process only one Data_type. Useful to rerun only nucleotide HGTs after query-ID matching fixes.')
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, 'jobs'), exist_ok=True)
    targets = pd.read_csv(args.targets, sep='\t', dtype=str)
    prot_tax = read_taxonomy_file(args.protein_taxonomy)
    nt_tax = read_taxonomy_file(args.nt_taxonomy)
    genome_taxid_map = read_genome_taxid_map(args.genome_taxid_map)

    print(f'[OK] loaded genome taxid map: {len(genome_taxid_map.get("by_genome", {}))} genomes')

    cache = {}
    for genome, dtype in targets[['Genome','Data_type']].drop_duplicates().itertuples(index=False):
        align_dir = args.protein_align_dir if dtype == 'Protein' else args.nt_align_dir
        out_file, tax_file = align_paths(align_dir, genome)
        key = (genome, dtype)
        try:
            cache[key] = load_hit_table(out_file, tax_file)
            print(f'[OK] loaded {dtype} hits for {genome}: {len(cache[key])} rows')
        except Exception as e:
            print(f'[WARN] cannot load {dtype} hits for {genome}: {e}')
            cache[key] = None

    manifest_rows = []
    all_selected_paths = []
    for _, row in targets.iterrows():
        hgt_id = row['hgt_id']
        dtype = row['Data_type']
        genome = row['Genome']
        jobdir = os.path.join(args.outdir, 'jobs', hgt_id)
        os.makedirs(jobdir, exist_ok=True)
        hits = cache.get((genome, dtype))
        if hits is None:
            meta = {'hgt_id':hgt_id, 'status':'missing_hit_file', 'selected_n':0, 'Data_type':dtype, 'Genome':genome, 'Accession':row.get('Accession',''), 'tree_priority':row.get('tree_priority','')}
            manifest_rows.append(meta)
            pd.DataFrame().to_csv(os.path.join(jobdir, 'selected_hits.tsv'), sep='\t', index=False)
            continue
        tax = prot_tax if dtype == 'Protein' else nt_tax
        selected, meta = sample_for_target(row, hits, tax, genome_taxid_map, args)
        out_path = os.path.join(jobdir, 'selected_hits.tsv')
        selected.to_csv(out_path, sep='\t', index=False)
        all_selected_paths.append(out_path)
        manifest_rows.append(meta)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(os.path.join(args.outdir, 'sampling_manifest.tsv'), sep='\t', index=False)

    nonempty = []
    for p in all_selected_paths:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                d = pd.read_csv(p, sep='\t', dtype=str)
                if not d.empty:
                    nonempty.append(d)
            except Exception:
                pass
    if nonempty:
        merged = pd.concat(nonempty, ignore_index=True)
        merged.to_csv(os.path.join(args.outdir, 'all_selected_hits.tsv'), sep='\t', index=False)
    else:
        pd.DataFrame().to_csv(os.path.join(args.outdir, 'all_selected_hits.tsv'), sep='\t', index=False)

    print(f'[OK] wrote sampling manifest to {os.path.join(args.outdir, "sampling_manifest.tsv")}')
    if not manifest.empty:
        print(manifest['status'].value_counts(dropna=False).to_string())
        cols = [c for c in ['self_n','close_n','out_host_n','out_other_n','unknown_n','selected_n'] if c in manifest.columns]
        if cols:
            print('\n[GROUP SUMMARY]')
            print(manifest[cols].apply(pd.to_numeric, errors='coerce').sum().to_string())


if __name__ == '__main__':
    main()
