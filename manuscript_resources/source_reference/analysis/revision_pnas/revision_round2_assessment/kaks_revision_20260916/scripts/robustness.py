from pathlib import Path
import json,re,hashlib
import pandas as pd
import numpy as np
from scipy.stats import spearmanr,rankdata
B=Path(__file__).resolve().parents[1]; R=B.parent; ROOT=R.parents[2]
T=B/'tables';rng=np.random.default_rng(20260916)
x=ROOT/'analysis/revision_pnas/article_before/Tables.xlsx'
s=pd.read_excel(x,sheet_name='Table S6',header=1).dropna(subset=['HGT_Accession'])
a=pd.read_excel(x,sheet_name='Table S5',header=1).dropna(subset=['Gene_ID'])
l=pd.read_csv(R/'innovation_final8989_20260915/inputs/final8989_record_locus_family_ledger.tsv',sep='\t')
m=pd.read_csv(R/'innovation_function_contribution_20260916/tables/species_metadata.tsv',sep='\t')
assert len(l)==8989 and l.accession.is_unique and s.HGT_Accession.is_unique and a.Gene_ID.is_unique
s=s.rename(columns={'HGT_Accession':'accession','Ka/Ks':'omega','Length.y':'aligned_nt'})
d=l.merge(s,on='accession',how='left',validate='one_to_one').merge(m[['species','clade','genus']],on='species',how='left',validate='many_to_one')
d=d.merge(a[['Gene_ID','Age_MYA','Species_Count']],left_on='accession',right_on='Gene_ID',how='left',validate='one_to_one')
d['age']=pd.to_numeric(d.original_minimal_age_Ma,errors='coerce')
for c in ['Ka','Ks','omega','aligned_nt']:d[c]=pd.to_numeric(d[c],errors='coerce')
d['strict_OG']=d.strict_own_OGs.fillna('').map(lambda z: z if re.fullmatch(r'OG\d+',z) else '')
d.loc[d.direct_evidence_status.ne('own_sequence_exact_eligible_representative'),'strict_OG']=''
d['old_S6_exact_match']=d.accession.isin(s.accession)
d['kaks_status']=np.where(d.old_S6_exact_match,'old_S6_saved_numeric','not_in_old_S6')
# New viral pairs are kept separate and never substituted for old pair partners.
p=pd.read_csv(ROOT/'analysis/revision_pnas/virus_supplement/downstream_novel/KaKs_pair_results.tsv',sep='\t')
n=p[p.accession.isin(l.loc[l.scope=='novel_viral_associated','accession']) & p.pair_role.eq('viral')]
assert len(n)==4
for _,q in n.iterrows():
 idx=d.index[d.accession.eq(q.accession)]; assert len(idx)==1 and not d.loc[idx,'old_S6_exact_match'].any()
 for src,dst in [('Ka','Ka'),('Ks','Ks'),('Ka/Ks','omega'),('Length','aligned_nt')]: d.loc[idx,dst]=q[src]
 d.loc[idx,'kaks_status']='new4_saved_viral_pair'
d['finite_positive_Ks']=np.isfinite(d[['Ka','Ks','omega']]).all(axis=1)&d.Ks.gt(0)&d.Ka.ge(0)&d.omega.ge(0)&d.aligned_nt.ge(150)
d['high_Ks_ge2']=d.Ks.ge(2)
# No full-sequence, orthology, donor or tree status is inferred here.
d.to_csv(T/'final8989_saved_estimate_status.tsv',sep='\t',index=False)
missing=s[~s.accession.isin(l.accession)].copy()
base=l.assign(base=l.accession.str.replace(r'\.\d+$','',regex=True)).groupby('base').accession.agg(lambda z:';'.join(z)).to_dict()
missing['same_base_catalogue_accessions']=missing.accession.str.replace(r'\.\d+$','',regex=True).map(base).fillna('')
missing['status']='unresolved_exact_accession_not_in_final_catalogue'
missing.to_csv(T/'unmatched_S6_accessions.tsv',sep='\t',index=False)
comp=d[d.Age_MYA.notna()&d.age.notna()].copy();comp['same_numeric_age']=np.isclose(comp.age,comp.Age_MYA,rtol=0,atol=1e-5)
comp[['record_id','age','Age_MYA','same_numeric_age']].to_csv(T/'age_field_concordance.tsv',sep='\t',index=False)

out=[]
def corr(z,x='age',y='omega'):
 if len(z)<3 or z[x].nunique()<2 or z[y].nunique()<2:return np.nan,np.nan
 c=spearmanr(z[x],z[y]);return float(c.statistic),float(c.pvalue)
def add(z,panel,unit,scope,nominal=False,boot=False):
 rho,pv=corr(z);row=dict(panel=panel,unit=unit,scope=scope,n=len(z),age_values=z.age.nunique(),rho=rho,naive_p_not_phylogenetic=(pv if nominal else np.nan),ci025=np.nan,ci975=np.nan)
 if boot and len(z)>=10 and np.isfinite(rho):
  arr=z[['age','omega']].to_numpy();b=[]
  for _ in range(1000):
   q=arr[rng.integers(0,len(arr),len(arr))]
   if len(np.unique(q[:,0]))>1 and len(np.unique(q[:,1]))>1:b.append(spearmanr(q[:,0],q[:,1]).statistic)
  row['ci025'],row['ci975']=np.quantile(b,[.025,.975]);row['valid_bootstrap_replicates']=len(b)
 out.append(row)
orig=s.merge(a[['Gene_ID','Age_MYA']],left_on='accession',right_on='Gene_ID',validate='one_to_one').rename(columns={'Age_MYA':'age'})
orig=orig[np.isfinite(orig[['age','omega']]).all(axis=1)]
add(orig,'submitted_S5_S6_including_zero','record','all',nominal=True)
add(orig[orig.age>0],'submitted_S5_S6_positive_age','record','all',nominal=True)
z=d[d.finite_positive_Ks&d.age.notna()]
add(z,'final_including_zero','record','all',nominal=True)
coverage=[];diagnostics=[]
for bound,label in [(None,'all_Ks'),(3,'Ks_lt3'),(2,'Ks_lt2'),(1,'Ks_lt1')]:
 q=z[z.age>0].copy()
 if bound is not None:q=q[q.Ks<bound]
 add(q,label,'record','all',nominal=True)
 coverage.append(dict(panel=label,records=len(q),species=q.species.nunique(),genera=q.genus.nunique(),clades=q.clade.nunique(),strict_records=q.strict_OG.ne('').sum(),global_strict_OGs=q.loc[q.strict_OG.ne(''),'strict_OG'].nunique()))
 for cl,sub in q.groupby('clade'):
  if len(sub)>=10:add(sub,label,'record',cl)
  left=q[q.clade!=cl]
  if len(left)>=10:add(left,label,'record','leave_out_'+cl)
 for y in ['Ka','Ks','omega']:
  rho,_=corr(q,y=y);diagnostics.append(dict(panel=label,comparison='age_vs_'+y,n=len(q),rho=rho))
 # Residual rank association removes species mean rank, not a phylogenetic model.
 rr=[]
 for sp,sub in q.groupby('species'):
  if len(sub)>=3 and sub.age.nunique()>=2 and sub.omega.nunique()>=2:
   ar=rankdata(sub.age);wr=rankdata(sub.omega)
   rr.extend(zip(ar-ar.mean(),wr-wr.mean()))
 if rr:
  r=np.asarray(rr);diagnostics.append(dict(panel=label,comparison='within_species_centered_rank_pearson',n=len(r),rho=np.corrcoef(r.T)[0,1]))
 t=q[q.strict_OG.ne('')]
 so=t.groupby(['species','clade','strict_OG'],as_index=False)[['age','omega']].median()
 add(so,label,'strict_species_OG_medians','all')
 go=so.groupby('strict_OG',as_index=False)[['age','omega']].median()
 add(go,label,'strict_global_OG_medians','all',boot=True)
 for cl,sub in so.groupby('clade'):
  gg=sub.groupby('strict_OG',as_index=False)[['age','omega']].median()
  if len(gg)>=10:add(gg,label,'strict_global_OG_medians',cl)
 for cl in sorted(so.clade.unique()):
  gg=so[so.clade!=cl].groupby('strict_OG',as_index=False)[['age','omega']].median()
  if len(gg)>=10:add(gg,label,'strict_global_OG_medians','leave_out_'+cl)
pd.DataFrame(out).to_csv(T/'correlation_sensitivity.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(T/'sensitivity_coverage.tsv',sep='\t',index=False)
pd.DataFrame(diagnostics).to_csv(T/'rate_component_diagnostics.tsv',sep='\t',index=False)
status=d.groupby(['scope','kaks_status','original_minimal_age_status'],dropna=False).size().rename('records').reset_index();status.to_csv(T/'catalogue_status_counts.tsv',sep='\t',index=False)
summary=dict(final_records=len(d),old_numeric=int(d.old_S6_exact_match.sum()),new_viral_numeric=len(n),no_saved_estimate=int(d.Ks.isna().sum()),age_both_available=len(comp),age_numeric_disagreement=int((~comp.same_numeric_age).sum()),s6_unmatched=len(missing),catalogue_unchanged=True,new_sequence_analysis=False)
(B/'robustness_status.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary));print(pd.DataFrame(out).query("scope=='all'").to_string(index=False))
