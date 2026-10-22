#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import re
import subprocess
import shutil
from contextlib import contextmanager
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Iterable
import pandas as pd

from hgt_tree_utils import read_fasta, write_fasta, crop_aa, translate_query_region, safe_id


def read_targets(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str)


def read_all_selected(jobs_dir: str) -> pd.DataFrame:
    frames = []
    for hgt_id in os.listdir(jobs_dir):
        p = os.path.join(jobs_dir, hgt_id, 'selected_hits.tsv')
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                frames.append(pd.read_csv(p, sep='\t', dtype=str))
            except Exception:
                pass
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    for c in ['qstart','qend','sstart','send','bitscore','aln_len','pident']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def recipient_self_taxid_from_hits(sub: pd.DataFrame) -> str:
    for col in ['recipient_self_taxid', 'self_taxid', 'hit_qtaxid']:
        if col in sub.columns and sub[col].notna().any():
            return str(sub[col].dropna().iloc[0])
    return 'NA'


def load_optional_fasta(path: str) -> Dict[str, str]:
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        return read_fasta(path)
    return {}


def write_id_file(ids: Iterable[str], path: str) -> int:
    ids = sorted(set(str(x) for x in ids if str(x) and str(x).lower() != 'nan'))
    with open(path, 'w') as out:
        for x in ids:
            out.write(x + '\n')
    return len(ids)


def samtools_faidx(genome_fasta: str, region: str) -> str:
    m = re.match(r'([^:]+):(\d+)-(\d+)$', region)
    if not m:
        raise ValueError(f'Cannot parse region: {region}')
    contig, a, b = m.group(1), int(m.group(2)), int(m.group(3))
    start, end = sorted([a, b])
    cmd = ['samtools', 'faidx', genome_fasta, f'{contig}:{start}-{end}']
    res = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return ''.join(line.strip() for line in res.stdout.splitlines() if not line.startswith('>'))


def parse_region_id(x: str) -> Optional[Tuple[str, int, int]]:
    m = re.match(r'^([^:]+):(\d+)-(\d+)$', str(x))
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def build_coord_index(fasta_dict: Dict[str, str]) -> Dict[str, List[Tuple[int, int, str]]]:
    idx: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for key in fasta_dict.keys():
        p = parse_region_id(key)
        if p:
            contig, a, b = p
            idx[contig].append((a, b, key))
    return idx


def get_nt_query_sequence(query_id: str, genome: str, nt_query_fasta: Dict[str,str], nt_coord_index: Dict[str, List[Tuple[int,int,str]]], genome_fasta_dir: str, coord_tolerance: int = 2) -> Optional[str]:
    # Exact first, because query IDs should usually match the curated Accession.
    if query_id in nt_query_fasta:
        return nt_query_fasta[query_id]
    for key in [query_id.replace(':','_'), query_id.replace(':','|')]:
        if key in nt_query_fasta:
            return nt_query_fasta[key]
    # Coordinate tolerant lookup for 0/1-based off-by-one differences, e.g. 175849 vs 175848.
    p = parse_region_id(query_id)
    if p:
        contig, a, b = p
        for aa, bb, key in nt_coord_index.get(contig, []):
            if abs(aa - a) <= coord_tolerance and abs(bb - b) <= coord_tolerance:
                return nt_query_fasta[key]
        # fallback to very high overlap if exact endpoint tolerance fails
        qlo, qhi = sorted([a, b])
        qlen = max(1, qhi - qlo + 1)
        best = None
        best_ov = 0.0
        for aa, bb, key in nt_coord_index.get(contig, []):
            rlo, rhi = sorted([aa, bb])
            ov = max(0, min(qhi, rhi) - max(qlo, rlo) + 1) / qlen
            if ov > best_ov:
                best_ov, best = ov, key
        if best is not None and best_ov >= 0.95:
            return nt_query_fasta[best]
    # Last fallback: extract from genome FASTA if provided.
    if genome_fasta_dir:
        genome_fasta = os.path.join(genome_fasta_dir, genome)
        if os.path.exists(genome_fasta):
            return samtools_faidx(genome_fasta, query_id)
    return None


@contextmanager
def open_fasta_stream(path: str, threads: int = 4):
    """Open plain or gzipped FASTA as text stream.

    For .gz, prefer pigz -dc -p threads when available; otherwise gzip -dc.
    This streams data and never loads the database FASTA into memory.
    """
    if str(path).endswith('.gz'):
        if shutil.which('pigz'):
            cmd = ['pigz', '-dc', '-p', str(max(1, int(threads))), path]
        else:
            cmd = ['gzip', '-dc', path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='replace')
        try:
            assert proc.stdout is not None
            yield proc.stdout
        finally:
            if proc.stdout:
                proc.stdout.close()
            err = proc.stderr.read() if proc.stderr else ''
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f'Command failed: {" ".join(cmd)}\n{err}')
    else:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            yield fh


def clean_aa(seq: str) -> str:
    return str(seq).replace('*','').replace('-','').replace('.', '').upper()


def make_subject_record_maps(selected: pd.DataFrame, args):
    """Return protein_map and nt_map: sseqid -> list of crop/write tasks.

    Unlike earlier versions, this function stores only crop instructions. The large database
    FASTA is streamed once and each matched subject is cropped immediately according to
    sstart/send before being written to a small cropped-fragment FASTA.
    """
    prot_map: Dict[str, List[dict]] = defaultdict(list)
    nt_map: Dict[str, List[dict]] = defaultdict(list)
    for _, r in selected.iterrows():
        sid = str(r.get('sseqid', ''))
        if not sid or sid.lower() == 'nan':
            continue
        dtype = str(r.get('data_type', r.get('Data_type', '')))
        task = {
            'hgt_id': str(r.get('hgt_id', '')),
            'sid': sid,
            'sstart': r.get('sstart'),
            'send': r.get('send'),
            'tree_role': str(r.get('tree_role', r.get('hit_group', 'HIT'))),
            'target_taxid': str(r.get('target_taxid', 'NA')),
            'species': str(r.get('species', 'species')),
        }
        if dtype == 'Protein':
            prot_map[sid].append(task)
        elif dtype == 'Nucleotide':
            # IMPORTANT: nucleotide translated DB IDs include frame suffix, e.g. XM_..._frame_1.
            # They must be matched exactly against FASTA first token to avoid pulling all six frames.
            nt_map[sid].append(task)
    return prot_map, nt_map


def stream_crop_subjects(fasta_path: str, records_by_id: Dict[str, List[dict]], out_fasta: str, summary_path: str, threads: int, flank_aa: int, label_prefix: str) -> None:
    if not fasta_path or not os.path.exists(fasta_path):
        raise FileNotFoundError(f'FASTA database not found: {fasta_path}')
    wanted = set(records_by_id.keys())
    if not wanted:
        open(out_fasta, 'w').close()
        pd.DataFrame([{'fasta': fasta_path, 'wanted_ids': 0, 'found_ids': 0, 'written_fragments': 0}]).to_csv(summary_path, sep='\t', index=False)
        return

    found_ids = set()
    written = 0
    cur_id = None
    keep = False
    chunks: List[str] = []

    def flush_current(out_handle):
        nonlocal written, chunks, cur_id, keep
        if not keep or cur_id is None:
            chunks = []
            return
        seq = ''.join(chunks)
        tasks = records_by_id.get(cur_id, [])
        if tasks:
            found_ids.add(cur_id)
        for task in tasks:
            frag = crop_aa(seq, task.get('sstart'), task.get('send'), flank_aa)
            frag = clean_aa(frag)
            if not frag:
                continue
            hgt_id = task['hgt_id']
            sid = task['sid']
            label = f"{task.get('tree_role','HIT')}|{hgt_id}|{safe_id(sid)}|taxid={safe_id(task.get('target_taxid','NA'))}|{safe_id(task.get('species','species'))}"
            out_handle.write(f'>{label}\n')
            for i in range(0, len(frag), 60):
                out_handle.write(frag[i:i+60] + '\n')
            written += 1
        chunks = []

    with open_fasta_stream(fasta_path, threads=threads) as fh, open(out_fasta, 'w') as out:
        for line in fh:
            if line.startswith('>'):
                flush_current(out)
                cur_id = line[1:].split()[0]
                # Exact first-token matching only. This is essential for *_frame_N IDs.
                keep = cur_id in wanted
                chunks = []
            else:
                if keep:
                    chunks.append(line.strip())
        flush_current(out)

    pd.DataFrame([{
        'fasta': fasta_path,
        'wanted_ids': len(wanted),
        'found_ids': len(found_ids),
        'missing_ids': len(wanted) - len(found_ids),
        'written_fragments': written,
        'out_fasta': out_fasta,
        'label_prefix': label_prefix,
        'match_mode': 'exact_first_token'
    }]).to_csv(summary_path, sep='\t', index=False)
    if len(wanted) - len(found_ids) > 0:
        miss_path = out_fasta + '.missing_ids.txt'
        with open(miss_path, 'w') as m:
            for sid in sorted(wanted - found_ids):
                m.write(sid + '\n')


def read_fragments_by_hgt(paths: List[str]) -> Dict[str, List[Tuple[str, str]]]:
    by_hgt: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    for path in paths:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        name = None
        chunks: List[str] = []
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if name is not None:
                        parts = name.split('|')
                        if len(parts) > 1:
                            by_hgt[parts[1]].append((name, ''.join(chunks)))
                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            if name is not None:
                parts = name.split('|')
                if len(parts) > 1:
                    by_hgt[parts[1]].append((name, ''.join(chunks)))
    return by_hgt


def best_query_window(hit_df: pd.DataFrame) -> Tuple[Optional[int], Optional[int]]:
    if hit_df.empty:
        return None, None
    vals = []
    for _, r in hit_df.iterrows():
        try:
            vals.append(int(r['qstart'])); vals.append(int(r['qend']))
        except Exception:
            pass
    if not vals:
        return None, None
    return min(vals), max(vals)


def build_one_final_fasta(hgt_row: pd.Series, hits: pd.DataFrame, fragments: Dict[str, List[Tuple[str,str]]], prot_query_fasta: Dict[str,str], nt_query_fasta: Dict[str,str], nt_coord_index: Dict[str, List[Tuple[int,int,str]]], args) -> Tuple[int, str]:
    hgt_id = str(hgt_row['hgt_id'])
    dtype = str(hgt_row['Data_type'])
    query_id = str(hgt_row['Accession'])
    genome = str(hgt_row['Genome'])
    jobdir = os.path.join(args.jobs_dir, hgt_id)
    os.makedirs(jobdir, exist_ok=True)
    sub = hits[hits['hgt_id'] == hgt_id].copy()
    if sub.empty:
        return 0, 'no_selected_hits'
    recipient_taxid = recipient_self_taxid_from_hits(sub)

    records: List[Tuple[str,str]] = []
    if dtype == 'Protein':
        qseq = prot_query_fasta.get(query_id)
        if qseq:
            qs, qe = best_query_window(sub)
            qseq2 = clean_aa(crop_aa(qseq, qs, qe, args.flank_aa))
            records.append((f'QUERY|{hgt_id}|{safe_id(query_id)}|taxid={safe_id(recipient_taxid)}', qseq2))
        else:
            with open(os.path.join(jobdir, 'query_missing.warning'), 'w') as w:
                w.write(f'Protein query sequence not found for {query_id}. Check --protein-query-fasta.\n')
    else:
        ntseq = get_nt_query_sequence(query_id, genome, nt_query_fasta, nt_coord_index, args.genome_fasta_dir, args.query_coord_tolerance)
        if ntseq:
            qs, qe = best_query_window(sub)
            if qs is not None:
                qaa = clean_aa(translate_query_region(ntseq, int(qs), int(qe), flank_nt=args.flank_aa*3))
                records.append((f'QUERY|{hgt_id}|{safe_id(query_id)}|taxid={safe_id(recipient_taxid)}', qaa))
        else:
            with open(os.path.join(jobdir, 'query_missing.warning'), 'w') as w:
                w.write(f'Nucleotide query sequence not found for {query_id}. Check --nt-query-fasta or --genome-fasta-dir.\n')

    records.extend(fragments.get(hgt_id, []))

    # Remove exact duplicate sequences and very short sequences.
    seen = set()
    clean_records = []
    for name, seq in records:
        seq = clean_aa(seq)
        if len(seq) < args.min_seq_len:
            continue
        if seq in seen:
            continue
        seen.add(seq)
        clean_records.append((name, seq))

    out_fa = os.path.join(jobdir, 'input.aa.fa')
    write_fasta(clean_records, out_fa)
    return len(clean_records), 'ok' if len(clean_records) >= 4 else 'too_few_sequences'


def make_cropped_subject_fastas(args, selected: pd.DataFrame):
    cache_dir = os.path.join(args.outdir, 'seq_cache')
    os.makedirs(cache_dir, exist_ok=True)
    prot_map, nt_map = make_subject_record_maps(selected, args)

    prot_ids = os.path.join(cache_dir, 'protein_subject_ids.exact.txt')
    nt_ids = os.path.join(cache_dir, 'nt_translated_subject_ids.exact.txt')
    n_prot = write_id_file(prot_map.keys(), prot_ids)
    n_nt = write_id_file(nt_map.keys(), nt_ids)

    prot_frag = os.path.join(cache_dir, 'protein_subject_fragments.cropped.fa')
    nt_frag = os.path.join(cache_dir, 'nt_translated_subject_fragments.cropped.fa')
    prot_summary = os.path.join(cache_dir, 'protein_subject_extraction_summary.tsv')
    nt_summary = os.path.join(cache_dir, 'nt_translated_subject_extraction_summary.tsv')

    if args.extract_sequences:
        print(f'[INFO] stream-crop protein subject IDs exactly: {n_prot}')
        stream_crop_subjects(args.protein_db_fasta, prot_map, prot_frag, prot_summary, args.decompress_threads, args.flank_aa, 'Protein')
        print(f'[INFO] stream-crop nt-translated subject IDs exactly: {n_nt}')
        stream_crop_subjects(args.nt_db_fasta, nt_map, nt_frag, nt_summary, args.decompress_threads, args.flank_aa, 'Nucleotide')
    else:
        print('[INFO] --extract-sequences not set; using existing cropped subject fragment FASTAs if present')
    return prot_frag, nt_frag


def main():
    ap = argparse.ArgumentParser(description='Extract/crop AA sequences for HGT tree FASTAs without writing huge full-sequence caches.')
    ap.add_argument('--targets', required=True)
    ap.add_argument('--jobs-dir', required=True)
    ap.add_argument('--protein-db', required=True, help='DIAMOND db path, kept for metadata/backward compatibility')
    ap.add_argument('--nt-db', required=True, help='DIAMOND db path, kept for metadata/backward compatibility')
    ap.add_argument('--protein-db-fasta', required=True, help='Original protein FASTA used to build the protein DIAMOND db')
    ap.add_argument('--nt-db-fasta', required=True, help='Original translated FASTA used to build the nt/protein DIAMOND db')
    ap.add_argument('--protein-query-fasta', default='')
    ap.add_argument('--nt-query-fasta', default='')
    ap.add_argument('--genome-fasta-dir', default='')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--flank-aa', type=int, default=20)
    ap.add_argument('--min-seq-len', type=int, default=40)
    ap.add_argument('--decompress-threads', type=int, default=8, help='Threads for pigz decompression of .gz FASTA if pigz is available')
    ap.add_argument('--query-coord-tolerance', type=int, default=2)
    ap.add_argument('--extract-sequences', action='store_true', help='Stream FASTA DBs once and write cropped subject fragments')
    # Backward compatibility: accepted but unused.
    ap.add_argument('--extractor', default='direct_crop')
    ap.add_argument('--seqkit-threads', type=int, default=4)
    ap.add_argument('--seqkit-match-name', type=int, default=0)
    ap.add_argument('--diamond-bin', default='diamond')
    ap.add_argument('--run-diamond', action='store_true', help='Deprecated alias for --extract-sequences')
    args = ap.parse_args()
    if args.run_diamond:
        args.extract_sequences = True

    targets = read_targets(args.targets)
    selected = read_all_selected(args.jobs_dir)
    if selected.empty:
        raise SystemExit('No selected_hits.tsv files found or all are empty.')

    prot_frag, nt_frag = make_cropped_subject_fastas(args, selected)
    fragments = read_fragments_by_hgt([prot_frag, nt_frag])

    prot_query_cache = load_optional_fasta(args.protein_query_fasta)
    nt_query_cache = load_optional_fasta(args.nt_query_fasta)
    nt_coord_index = build_coord_index(nt_query_cache)

    manifest = []
    for _, row in targets.iterrows():
        hgt_id = str(row['hgt_id'])
        n, status = build_one_final_fasta(row, selected, fragments, prot_query_cache, nt_query_cache, nt_coord_index, args)
        manifest.append({
            'hgt_id': hgt_id,
            'Data_type': row['Data_type'],
            'Accession': row['Accession'],
            'tree_priority': row.get('tree_priority',''),
            'seq_n': n,
            'status': status
        })
    pd.DataFrame(manifest).to_csv(os.path.join(args.outdir, 'fasta_manifest.tsv'), sep='\t', index=False)
    print(f'[OK] wrote cropped FASTAs under {args.jobs_dir}; manifest: {os.path.join(args.outdir, "fasta_manifest.tsv")}')


if __name__ == '__main__':
    main()
