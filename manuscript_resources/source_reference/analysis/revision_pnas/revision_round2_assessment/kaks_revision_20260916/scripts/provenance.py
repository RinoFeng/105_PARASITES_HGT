from pathlib import Path
import pandas as pd,json,hashlib
B=Path(__file__).resolve().parents[1];R=B.parent;ROOT=R.parents[2];T=B/'tables'
x=ROOT/'analysis/revision_pnas/article_before/Tables.xlsx'
s=pd.read_excel(x,sheet_name='Table S6',header=1).dropna(subset=['HGT_Accession'])
l=pd.read_csv(R/'innovation_final8989_20260915/inputs/final8989_record_locus_family_ledger.tsv',sep='\t')
legacy=ROOT/'analysis/Results_article_new/KaKs_KnKs/prot/all_HGT.prot.kaks'
idx={}
with legacy.open() as f:
 for line in f:
  p=line.rstrip('\n').split('\t')
  if len(p)<5 or p[1]=='Method':continue
  key=p[0].split('/',1)[-1]
  try:v=tuple(float(t) for t in p[2:5])
  except ValueError:continue
  idx.setdefault(key,[]).append((p[1],v))
rows=[]
import numpy as np
for _,q in s.iterrows():
 matches=idx.get(q.Sequence_pair,[])
 good=[mo for mo,v in matches if np.allclose(v,[q.Ka,q.Ks,q['Ka/Ks']],rtol=1e-7,atol=1e-9)]
 rows.append(dict(accession=q.HGT_Accession,exact_pair_matches=len(matches),numeric_matching_models=';'.join(sorted(set(good))),status='matched' if good else 'unresolved'))
d=pd.DataFrame(rows);d.to_csv(T/'S6_legacy_model_provenance.tsv',sep='\t',index=False)
# Investigate identifiers without substituting a different candidate's estimate.
a=pd.read_csv(R/'innovation_function_contribution_20260916/cache/accession_index.tsv.gz',sep='\t')
s2=pd.read_excel(x,sheet_name='Table S2',header=1)
u=pd.read_csv(T/'unmatched_S6_accessions.tsv',sep='\t')
u['present_in_submitted_S2']=u.accession.isin(s2.HGT_Accession)
u['present_in_current_proteome_index']=u.accession.isin(a.accession)
u['status']='author_confirmed_intentional_PCR_candidate_replacement; do_not_transfer_values'
u.to_csv(T/'unmatched_S6_accessions.tsv',sep='\t',index=False)
a[a.accession.isin(u.accession)].to_csv(T/'unmatched_S6_current_input_links.tsv',sep='\t',index=False)
old_catalogue=l[l.scope=='original_8985']
so=s2[~s2.HGT_Accession.isin(old_catalogue.accession)]
lo=old_catalogue[~old_catalogue.accession.isin(s2.HGT_Accession)]
so.to_csv(T/'submitted_S2_not_exact_current_old.tsv',sep='\t',index=False)
lo.to_csv(T/'current_old_not_exact_submitted_S2.tsv',sep='\t',index=False)
fc=pd.read_csv(R/'innovation_final8989_20260915/inputs/final_catalogue_8989.tsv',sep='\t',usecols=['Specie','Accession','Position(Start-End)'])
fc=fc[fc.Accession.isin(lo.accession)]
cc=so[['Specie','HGT_Accession','Position(Start-End)']].merge(fc,on='Specie',validate='one_to_one',suffixes=('_submitted','_current'))
import re
def coords(v):return tuple(sorted(int(z) for z in re.findall(r'\d+',str(v))))
cc['same_coordinate_numbers']=cc['Position(Start-End)_submitted'].map(coords)==cc['Position(Start-End)_current'].map(coords)
cc['interpretation']='same-species ID differences; not an approved alias mapping'
cc.to_csv(T/'four_ID_coordinate_comparison.tsv',sep='\t',index=False)

files=[x,legacy,R/'innovation_final8989_20260915/inputs/final_catalogue_8989.tsv',R/'innovation_final8989_20260915/inputs/final8989_record_locus_family_ledger.tsv',R/'innovation_function_contribution_20260916/cache/observations.tsv.gz',R/'innovation_execution/family_comparison/resolved_gene_units.tsv.gz']
manifest=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]
(B/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
summary=dict(submitted_rows=len(d),numeric_matched_to_saved_MA=int(d.numeric_matching_models.eq('MA').sum()),unresolved_model_provenance=int(d.status.ne('matched').sum()),old_software_executable_version='not established from saved results',new_viral_software_record='MUSCLE 3.8.31; saved KaKs model MA',unmatched_accessions=len(u),unmatched_in_submitted_S2=int(u.present_in_submitted_S2.sum()),unmatched_in_current_proteome=int(u.present_in_current_proteome_index.sum()),id_change_relationship='author-confirmed intentional restoration of PCR-tested candidates; current catalogue authoritative; old values not reassigned')
(B/'provenance_status.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
