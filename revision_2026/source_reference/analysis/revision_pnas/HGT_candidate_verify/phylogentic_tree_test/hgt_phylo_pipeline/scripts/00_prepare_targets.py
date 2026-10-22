#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import pandas as pd
import re


def safe_hgt_id(row, idx: int) -> str:
    sp = re.sub(r'[^A-Za-z0-9]+', '_', str(row.get('Specie', 'Species')))[:40]
    acc = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(row.get('Accession', 'ACC')))[:80]
    return f'HGT{idx:06d}__{sp}__{acc}'


def main():
    ap = argparse.ArgumentParser(description='Prepare curated HGT targets for tree construction.')
    ap.add_argument('--curated-tsv', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--keep-main', nargs='+', default=['Curated_main_set', 'Tier3_rescued_by_extra_support'])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.curated_tsv, sep='\t', dtype=str)
    df = df[df['main_analysis_set'].isin(args.keep_main)].copy()
    df = df.reset_index(drop=True)
    df['hgt_id'] = [safe_hgt_id(row, i+1) for i, row in df.iterrows()]

    # Normalize numeric fields used downstream.
    for col in ['Self_num','Close_num','Out_num','HGT_coverage','HGT_identity','alien_index','Length','Entropy','score']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    keep_cols = [
        'hgt_id','Specie','Genome','Accession','query','Data_type','Target_taxid','Target_lineage',
        'Out_pct','Self_num','Close_num','Out_num','HGT_coverage','HGT_identity','alien_index','Length',
        'Entropy','main_analysis_set','contig_context_tier','phylo_testability','tree_priority',
        'target_domain','data_type_overlap','Parasitism_type','Parasite_clade','Description','Preferred_name'
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    out = df[keep_cols].copy()
    out.to_csv(os.path.join(args.outdir, 'targets_for_tree.tsv'), sep='\t', index=False)

    summary = out.groupby(['Data_type','tree_priority'], dropna=False).size().reset_index(name='n')
    summary.to_csv(os.path.join(args.outdir, 'targets_for_tree_summary.tsv'), sep='\t', index=False)
    print(f'[OK] Wrote {len(out)} targets to {os.path.join(args.outdir, "targets_for_tree.tsv")}')
    print(summary.to_string(index=False))

if __name__ == '__main__':
    main()
