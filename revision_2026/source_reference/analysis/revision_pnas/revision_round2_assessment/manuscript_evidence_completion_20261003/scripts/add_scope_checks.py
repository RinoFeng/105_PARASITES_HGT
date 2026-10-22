from pathlib import Path
import pandas as pd,numpy as np,json
O=Path(__file__).resolve().parents[1];D=O/'data';R=O.parent
x=pd.read_csv(D/'age_record_sensitivity.tsv',sep='\t');
def era(v):return 'unresolved_zero' if v<=0 else 'Cenozoic' if v<=66 else 'Mesozoic' if v<=252 else 'Paleozoic' if v<=541 else 'Precambrian'
# Match the original R boundaries (541 Ma), retaining original figure exclusions explicitly.
rows=[]
for col,label in [('raw_MRCA_Ma','raw_all_occurrences'),('saved_age_Ma','one_pass_LOO')]:
 for scope,z in [('all_saved_age_records',x),('positive_age',x[x[col]>0]),('original_plot_rule',x[(x[col]>0)&(x[col]<=1500)&((x[col]-1448.548).abs()>1e-6)]),('fixed_saved_plot_subset',x[(x.saved_age_Ma>0)&(x.saved_age_Ma<=1500)&((x.saved_age_Ma-1448.548).abs()>1e-6)])]:
  for e,n in z[col].map(era).value_counts().items():rows.append(dict(version=label,scope=scope,era=e,records=n,denominator=len(z),fraction=n/len(z),median_Ma=z[col].median()))
pd.DataFrame(rows).to_csv(D/'age_plot_scope_sensitivity.tsv',sep='\t',index=False)
x['raw_era']=x.raw_MRCA_Ma.map(era);x['saved_era']=x.saved_age_Ma.map(era)
x.groupby(['Data_type','raw_era','saved_era']).size().rename('records').reset_index().to_csv(D/'age_era_transitions.tsv',sep='\t',index=False)
rows=[]
for col,label in [('raw_MRCA_Ma','raw_all_occurrences'),('saved_age_Ma','one_pass_LOO')]:
 for typ,z in [('ALL',x)]+list(x.groupby('Data_type')):
  for e,n in z[col].map(era).value_counts().items():rows.append(dict(version=label,data_type=typ,era=e,records=n,denominator=len(z),fraction=n/len(z)))
pd.DataFrame(rows).to_csv(D/'age_distribution_before_after.tsv',sep='\t',index=False)
u=pd.read_csv(D/'final48_locus_inventory.tsv.gz',sep='\t');p=pd.read_excel(R/'manuscript_completed_20260929/tables/Table_S2_Family_abundance.xlsx',sheet_name='Q_panel48')
v=u.drop_duplicates(['species','contig','strand','coding_structure']);counts=v.groupby(['species','role']).agg(families=('OG','nunique'),loci=('unit_id','nunique')).reset_index();rows=[]
for r in p.itertuples():
 c=counts[(counts.species==r.species)&counts.role.eq('candidate')].iloc[0];b=counts[(counts.species==r.species)&counts.role.eq('background')].iloc[0]
 old=(r.N/r.C)/(r.M/r.B);new=(c.loci/c.families)/(b.loci/b.families)
 rows.append(dict(species=r.species,original_N=r.N,original_M=r.M,coordinate_N=c.loci,coordinate_M=b.loci,original_C=r.C,original_B=r.B,coordinate_C=c.families,coordinate_B=b.families,original_within_ratio=old,coordinate_within_ratio=new,relative_ratio_change=new/old-1))
pd.DataFrame(rows).to_csv(D/'coordinate_dedup_sensitivity.tsv',sep='\t',index=False)
print('positive',int((x.saved_age_Ma>0).sum()),'plot',(x.saved_age_Ma.gt(0)&(x.saved_age_Ma-1448.548).abs().gt(1e-6)).sum(),'saved_root',np.isclose(x.saved_age_Ma,1448.548,atol=1e-6,rtol=0).sum());print(pd.DataFrame(rows).query('relative_ratio_change != 0')[['species','original_M','coordinate_M','relative_ratio_change']].to_string(index=False))
