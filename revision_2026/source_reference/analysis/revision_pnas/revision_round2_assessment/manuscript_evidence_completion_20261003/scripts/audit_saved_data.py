"""Audit saved tables and numerical tree products; never read runtime logs."""
from pathlib import Path
from functools import lru_cache
import hashlib,json,re,csv
import pandas as pd
import numpy as np
from Bio import Phylo
ROOT=Path(__file__).resolve().parents[1]; PROJECT=ROOT.parents[3]; R=ROOT.parent
CURRENT=R/'manuscript_completed_20260929'; A=PROJECT/'analysis/Results_article_new/HGT_time_count';V=PROJECT/'analysis/revision_pnas/virus_supplement'
T=ROOT/'data'; manifest=[]
def register(p):
 manifest.append(dict(path=str(p.resolve()),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return p
def write(name,df):df.to_csv(T/name,sep='\t',index=False)
def readtab(file,sheet):return pd.read_excel(register(CURRENT/'tables'/file),sheet_name=sheet)
# Reproduce the actual one-pass leave-one-species-out procedure on saved occurrences.
tree=Phylo.read(register(A/'phylogenetic_tree/105_parasites_timetree.timetree.nex'),'nexus'); tips={t.name:t for t in tree.get_terminals()}
@lru_cache(None)
def age(names):
 if len(names)<2:return 0.
 node=tree.common_ancestor(list(names));return tree.distance(node,node.get_terminals()[0])
@lru_cache(None)
def loo(names):
 raw=age(names);best=raw;omitted=[]
 if raw>=200 and len(names)>2:
  vals={s:age(tuple(t for t in names if t!=s)) for s in names}
  best=min([raw]+[x for x in vals.values() if x<raw*.5]);omitted=[s for s,x in vals.items() if x==best and best<raw]
 return raw,best,';'.join(omitted)
x=readtab('Table_S6_Historical_depth.xlsx','Minimum_age_occurrence').fillna('');old=[]
for typ,fn in [('Protein','Final_HGT_Robust_TBLASTN.csv'),('Nucleotide','Final_HGT_Robust_TBLASTX.csv')]:
 z=pd.read_csv(register(A/fn)).fillna('');z['Data_type']=typ;old.append(z)
old=pd.concat(old,ignore_index=True);assert not old.duplicated(['Gene_ID','Data_type']).any()
z=x.merge(old,on=['Gene_ID','Data_type'],suffixes=('_table','_saved'),how='left',validate='many_to_one',indicator=True)
rows=[]
for r in z.to_dict('records'):
 names=tuple(sorted(set(str(r['All_Species_List_table']).split(';'))-{''}));missing=set(names)-tips.keys()
 raw,best,possible=(np.nan,np.nan,'') if missing else loo(names)
 saved=float(r['Age_MYA_table']);removed=str(r['Removed_Outliers_table']).strip();removed='' if removed.lower() in {'none','nan','na','-'} else removed
 rows.append(dict(Gene_ID=r['Gene_ID'],Data_type=r['Data_type'],saved_age_Ma=saved,raw_MRCA_Ma=raw,recomputed_LOO_age_Ma=best,changed=bool(np.isfinite(raw) and abs(raw-best)>1e-4),saved_reproduced=bool(np.isclose(saved,best,atol=.005,rtol=0)),source_row_matches=bool(r['_merge']=='both' and np.isclose(saved,float(r['Age_MYA_saved']),atol=1e-6)),species_count=len(names),recorded_species_count=r['Species_Count'],removed_species=removed,permitted_single_removals=possible,removed_label_consistent=(removed in possible.split(';') if possible else not removed),raw_members= ';'.join(names),remaining_members=';'.join(s for s in names if s!=removed),missing_tree_names=';'.join(sorted(missing)),removal_evidence='age_change_only_no_independent_rejection_in_script' if removed else 'not_applicable'))
y=pd.DataFrame(rows);write('age_record_sensitivity.tsv',y)
def era(v):
 return 'unresolved_zero' if v<=0 else 'Cenozoic' if v<66 else 'Mesozoic' if v<252 else 'Paleozoic' if v<=540 else 'Precambrian'
y['raw_era']=y.raw_MRCA_Ma.map(era);y['saved_era']=y.saved_age_Ma.map(era);write('age_era_transitions.tsv',y.groupby(['Data_type','raw_era','saved_era'],dropna=False).size().rename('records').reset_index())
dist=[]
for label,cols in [('raw_all_occurrences','raw_MRCA_Ma'),('one_pass_LOO','saved_age_Ma')]:
 for typ,zz in [('ALL',y)]+list(y.groupby('Data_type')):
  for e,n in zz[cols].map(era).value_counts().items():dist.append(dict(version=label,data_type=typ,era=e,records=n,denominator=len(zz),fraction=n/len(zz)))
write('age_distribution_before_after.tsv',pd.DataFrame(dist))
# Calibrations: compare saved input with recovered tree; do not infer a release provenance.
cals=[]
for line in register(A/'phylogenetic_tree/perfect_date.txt').read_text().splitlines():
 if not line.strip():continue
 pair,date=line.split();names=tuple(pair.split(','));ok=all(n in tips for n in names);raw=age(tuple(sorted(names))) if ok else np.nan
 cals.append(dict(taxon_a=names[0],taxon_b=names[1],input_date=float(date),input_age_Ma=-float(date),constraint='point_value_in_perfect_date.txt',both_tips_present=ok,tree_node_age_Ma=raw,difference_Ma=raw+float(date),source='TimeTree as recorded in original input; independent export not supplied'))
write('calibration_input_audit.tsv',pd.DataFrame(cals))
# Viral whole-catalogue coverage, not a new homology search.
v=readtab('Table_S3_Source_evidence.xlsx','Catalogue_viral_evidence').fillna('');v['reviewed']=v.uniform_search_reviewed.astype(str).str.lower().eq('true')
write('viral_ranking_coverage.tsv',v.groupby(['cohort_role','reviewed','uniform_viral_rank_status'],dropna=False).size().rename('records').reset_index())
write('viral_ranking_per_record.tsv',v[['hgt_id','Specie','Accession','cohort_role','reviewed','uniform_viral_rank_status','uniform_unknown_competes','raw_query_comparability','display_best_distant_domain','display_best_distant_basis']])
# Exact original plot inputs, verifying row order by original target taxonomy IDs.
p=register(PROJECT/'analysis/revision_pnas/HGT_candidate_verify/main_phylogeny_retained_HGT_8985.tsv');h=pd.read_csv(p,sep='\t',keep_default_na=False)
p=register(PROJECT/'analysis/revision_pnas/HGT_candidate_verify/HGT_target_lineage.txt');lineage=pd.read_csv(p,sep='\t',header=None,keep_default_na=False,dtype=str)
assert len(h)==len(lineage)==8985
print('original columns',h.columns.tolist(),flush=True)
ph=lineage.iloc[:,2].map(lambda s:s.split(';')[1].strip() if len(s.split(';'))>1 and s.split(';')[1].strip() else 'Unclassified');h['Phylum']=ph
h['Parasite_clade_display']=h.Parasite_clade.replace({'Protozoa':'Microbial eukaryotes','microbial eukaryotes':'Microbial eukaryotes','Metazoa':'Metazoan parasites'})
h['Parasitism_type']=h.Parasitism_type.map(lambda s:'Intracellular' if '胞内寄生' in s or 'Intracellular' in s else 'Extracellular' if '胞外寄生' in s or 'Extracellular' in s else 'Free-living' if '自由生活' in s or 'Free-living' in s else 'Unclassified')
h=h[h.Parasitism_type.ne('Unclassified')&h.Specie.ne('')].copy();h['Niche']=h.Parasitism_type+' '+h.Parasite_clade_display
assert h.groupby('Specie').Niche.nunique().max()==1
counts=pd.crosstab(h.Specie,h.Phylum).sort_index();write('dispersion_phylum_counts.tsv',counts.reset_index());write('dispersion_metadata.tsv',h[['Specie','Niche','Parasite_clade_display','Parasitism_type']].drop_duplicates().sort_values('Specie'))
# Exact identical-sequence units survive gene-level deduplication.
u=pd.read_csv(register(R/'innovation_execution/family_comparison/resolved_gene_units.tsv.gz'),sep='\t',usecols=['species','genome','gene_id','sequence_id','unit_id'],keep_default_na=False)
assert not u.duplicated(['genome','gene_id']).any();assert u.unit_id.is_unique
blocks=u.groupby(['genome','sequence_id']).agg(gene_units=('gene_id','nunique'),unit_ids=('unit_id',lambda x:';'.join(x))).reset_index();same=blocks[blocks.gene_units>1].copy();write('same_sequence_distinct_gene_units.tsv',same)
# Family panel quantities are reconstructed from the declared family-count table.
panel=readtab('Table_S2_Family_abundance.xlsx','Q_panel48');fam=pd.read_csv(register(R/'metazoa_extension_20260922/tables/species_family_counts.tsv.gz'),sep='\t');fam=fam[(fam.scenario=='augmented_same_reference')&fam.species.isin(panel.species)]
summary=fam.groupby(['species','role']).agg(families=('OG','nunique'),loci=('loci','sum')).reset_index();write('family_panel_count_reconstruction.tsv',summary)
checks=[]
for row in panel.itertuples():
 for role,fc,lc in [('candidate','C','N'),('background','B','M')]:
  rr=summary[(summary.species==row.species)&(summary.role==role)].iloc[0];checks.append(dict(species=row.species,role=role,family_match=int(rr.families)==getattr(row,fc),locus_match=int(rr.loci)==getattr(row,lc)))
write('family_panel_count_checks.tsv',pd.DataFrame(checks))
# Download snapshots: cryptographic values preserved verbatim, independent check performed below.
resources=json.loads(register(V/'databases/download_manifest.json').read_text());db=[]
for rec in resources:
 p=V/'databases'/rec['file'];sha=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):sha.update(b)
 db.append(dict(**rec,sha256_current=sha.hexdigest(),hash_matches=sha.hexdigest()==rec['sha256'],release_number='not_recorded'))
write('database_snapshots.tsv',pd.DataFrame(db))
result=dict(age_records=len(y),saved_ages_reproduced=int(y.saved_reproduced.sum()),source_rows_match=int(y.source_row_matches.sum()),changed_ages=int(y.changed.sum()),era_changed=int((y.raw_era!=y.saved_era).sum()),removal_labels_consistent=int(y.removed_label_consistent.sum()),raw_median=float(y.raw_MRCA_Ma.median()),saved_median=float(y.saved_age_Ma.median()),tree_tips=len(tips),tree_root_depth_min=min(tree.depths()[t] for t in tree.get_terminals()),tree_root_depth_max=max(tree.depths()[t] for t in tree.get_terminals()),calibration_input_rows=len(cals),dated_tree_rebuilt=False,sequence_search_rerun=False,viral_total=len(v),viral_reviewed=int(v.reviewed.sum()),viral_reviewed_status=v[v.reviewed].uniform_viral_rank_status.value_counts().to_dict(),gene_units=len(u),identical_sequence_multilocus_blocks=len(same),genes_in_identical_sequence_blocks=int(same.gene_units.sum()),family_panel_checks_pass=all(c['family_match'] and c['locus_match'] for c in checks),database_hashes_pass=all(d['hash_matches'] for d in db))
(ROOT/'checks/saved_data_summary.json').write_text(json.dumps(result,indent=2)+'\n');write('input_manifest.tsv',pd.DataFrame(manifest).drop_duplicates('path'));print(json.dumps(result,indent=2),flush=True)
