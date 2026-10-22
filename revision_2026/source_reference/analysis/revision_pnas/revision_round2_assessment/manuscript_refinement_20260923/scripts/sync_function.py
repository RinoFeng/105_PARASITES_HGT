from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[1];R=B.parent;M=R/'metazoa_extension_20260922';T=B/'tables'
CATS=list('EFGH');PARTS=CATS+['Other','Unknown'];CL=['Amoebozoa','Apicomplexa','Euglenozoa','Metazoa'];checks=[];manifest=[]
def read(p):
 manifest.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return pd.read_csv(p,sep='\t')
def save(x,n):x.to_csv(T/n,sep='\t',index=False)
def check(n,v):
 checks.append({'check':n,'passed':bool(v)});assert v,n
f=read(M/'phylogenetic_comparison/family_inputs48.tsv.gz').set_index(['species','OG','role'])
pairs=read(M/'tables/matched_pairs.tsv.gz');old_effect=read(M/'tables/species_function_effects.tsv');base_species=read(M/'tables/species_quantity.tsv')
allc=[];alld=[];alls=[];allcurves=[];allboot=[];allcl=[];allprimary=[];deldefs=[]
for sc in ['original47','augmented_same_reference']:
 pp=pairs[pairs.scenario==sc];rows=[]
 for (sp,o),g in pp.groupby(['species','candidate_OG'],sort=True):
  c=f.loc[(sp,o,'candidate')];bg=f.loc[[(sp,x,'background') for x in g.background_OG]]
  row={'species':sp,'clade':c.clade,'OG':o,'loci':c.units}
  for k in PARTS:row['candidate_'+k]=c[k];row['background_'+k]=bg[k].mean();row[k]=c[k]-bg[k].mean()
  rows.append(row)
 a=pd.DataFrame(rows);nsp=a.groupby('clade').species.nunique();nf=a.groupby('species').size()
 a['weight']=1/a.species.map(nf)/a.clade.map(nsp)/4
 check(sc+'_matched_units',len(a)==(1929 if sc=='original47' else 1940))
 check(sc+'_weight_sum',np.isclose(a.weight.sum(),1))
 ss=[]
 for (sp,cl),z in a.groupby(['species','clade']):
  for k in PARTS:ss.append({'species':sp,'clade':cl,'category':k,'candidate_share':z['candidate_'+k].mean(),'background_share':z['background_'+k].mean(),'effect':z[k].mean(),'matched_families':len(z),'scenario':sc})
 species=pd.DataFrame(ss); prior=old_effect[old_effect.scenario==sc]
 j=species.merge(prior,on=['species','category'],suffixes=('_new','_old'))
 check(sc+'_prior_effects_reproduced',np.allclose(j.effect_new,j.effect_old))
 allcl.append(species.groupby(['scenario','clade','category'],as_index=False)[['candidate_share','background_share','effect']].mean());alls.append(species)
 cr=[]
 for k in CATS:
  z=a[['species','clade','OG','weight','loci',k]].copy().rename(columns={k:'difference'});z['category']=k;z['scenario']=sc
  z['primary_contribution']=z.weight*z.difference;z['primary_positive_part']=z.weight*z.difference.clip(lower=0);z['primary_negative_part']=z.weight*z.difference.clip(upper=0);cr.append(z)
 con=pd.concat(cr,ignore_index=True);byog=con.groupby(['scenario','category','OG'],as_index=False)[['primary_contribution','primary_positive_part','primary_negative_part']].sum();allc.append(con)
 for k,z in byog.groupby('category'):
  v=z[z.primary_positive_part>0].sort_values(['primary_positive_part','OG'],ascending=[False,True]).copy();v['rank']=np.arange(1,len(v)+1);v['cumulative_positive_share']=v.primary_positive_part.cumsum()/v.primary_positive_part.sum();v['fraction_of_positive_OGs']=v['rank']/len(v);allcurves.append(v)
  shares=v.primary_positive_part/v.primary_positive_part.sum()
  summary={'scenario':sc,'category':k,'positive_OGs':len(v),'n50':int(v.loc[v.cumulative_positive_share>=.5,'rank'].iloc[0]),'n80':int(v.loc[v.cumulative_positive_share>=.8,'rank'].iloc[0]),'top5_share':shares.head(5).sum(),'gross_positive_pp':z.primary_positive_part.sum()*100,'gross_negative_pp':z.primary_negative_part.sum()*100,'net_effect_pp':z.primary_contribution.sum()*100}
  alld.append(summary)
  # Largest-locus and largest-positive-contribution rules use current matched panel.
  loc=a.groupby('OG').loci.sum().reset_index().sort_values(['loci','OG'],ascending=[False,True]); sets={'remove_top5_locus_OGs':set(loc.head(5).OG),'remove_top5_category_positive_OGs':set(v.head(5).OG)}
  for name,remove in sets.items():
   deldefs.append({'scenario':sc,'category':k,'deletion':name,'OGs':';'.join(sorted(remove))})
   rem=a[~a.OG.isin(remove)].groupby(['species','clade']).agg(n=('OG','size'),effect=(k,'mean')).reset_index();rem=rem[rem.n>=10];nc=rem.groupby('clade').species.nunique();clades=list(nc[nc>=3].index);rem=rem[rem.clade.isin(clades)]
   baseline=species[(species.category==k)&species.species.isin(rem.species)]
   deletion={'scenario':sc,'category':k,'deletion':name,'retained_species':rem.species.nunique(),'retained_clades':len(clades),'baseline_same_species_pp':baseline.groupby('clade').effect.mean().mean()*100,'remaining_effect_pp':rem.groupby('clade').effect.mean().mean()*100,'baseline_all_species_pp':summary['net_effect_pp'],'comparable':len(clades)>=3}
   allboot.append(deletion)
 # Original shared candidate-OG Poisson weights; controls remain fixed conditional expectations.
 rng=np.random.default_rng(20260923);ogs=sorted(a.OG.unique());oi={o:i for i,o in enumerate(ogs)};mult=rng.poisson(1,(1000,len(ogs)))[:,a.OG.map(oi).to_numpy()];bootsp={}
 for sp,z in a.groupby('species'):
  idx=z.index.to_numpy();mm=mult[:,idx];den=mm.sum(axis=1);vals=np.divide(mm@z[CATS].to_numpy(),den[:,None],out=np.full((1000,4),np.nan),where=den[:,None]>0);bootsp[sp]=vals
 bootcl={}
 for cl in CL:
  sps=sorted(a.loc[a.clade==cl,'species'].unique());bs=np.nanmean([bootsp[s] for s in sps],axis=0);bootcl[cl]=bs
  for ki,k in enumerate(CATS):
   ci=np.nanquantile(bs[:,ki],[.025,.975]);allprimary.append({'scenario':sc,'scope':cl,'category':k,'q025':ci[0],'q975':ci[1],'replicates':1000,'seed':20260923})
 bp=np.mean(list(bootcl.values()),axis=0)
 for ki,k in enumerate(CATS):
  ci=np.nanquantile(bp[:,ki],[.025,.975]);allprimary.append({'scenario':sc,'scope':'equal_four_clades','category':k,'q025':ci[0],'q975':ci[1],'replicates':1000,'seed':20260923})
 # Keep family rows for reproducibility, including Other/Unknown denominator.
 a['scenario']=sc;save(a,'family_effects_'+sc+'.tsv.gz')
 if sc=='original47':
  old=read(R/'family_contribution_extension_20260922/tables/contribution_summary.tsv').set_index('category')
  for row in [d for d in alld if d['scenario']==sc]:
   for k in ['positive_OGs','n50','top5_share','net_effect_pp']:check('original47_'+row['category']+'_'+k,np.isclose(row[k],old.loc[row['category'],k]))
save(pd.concat(allc),'family_contributions_47_48.tsv.gz');save(pd.DataFrame(alld),'contribution_summary_47_48.tsv');save(pd.concat(allcurves),'contribution_curves_47_48.tsv');save(pd.concat(alls),'species_functions_47_48.tsv');save(pd.concat(allcl),'clade_functions_47_48.tsv');save(pd.DataFrame(allboot),'deletion_same_scope_47_48.tsv');save(pd.DataFrame(deldefs),'deletion_definitions_47_48.tsv');save(pd.DataFrame(allprimary),'shared_OG_intervals_47_48.tsv')
primary=pd.concat(allcl).groupby(['scenario','category'],as_index=False)[['candidate_share','background_share','effect']].mean();save(primary,'primary_functions_47_48.tsv')
check('effect_sums',np.allclose(primary.groupby('scenario').effect.sum(),0,atol=1e-12))
# Raw components derived using the already recorded equal-size species subsets.
draw=read(M/'tables/equal_species_draws.tsv.gz');sq=base_species[base_species.scenario=='augmented_same_reference'].set_index('species');dr=draw[draw.scenario=='augmented_same_reference'].copy()
for idx,r in dr.iterrows():
 z=sq.loc[r.species.split(';')];dr.loc[idx,'raw_family']=z.C.mean();dr.loc[idx,'raw_within']=z.N.sum()/z.C.sum()
qq=[];ref=dr[dr.clade=='Metazoa'].set_index('replicate')
for cl in CL[:-1]:
 z=dr[dr.clade==cl].set_index('replicate')
 for version,cols in [('raw',['raw_family','raw_within']),('ordinary_normalized',['family_ratio','within_family_loci_ratio'])]:
  for col in cols:
   vals=np.log2(z[col]/ref[col]);qq.append({'clade':cl,'version':version,'component':'family' if col in ['raw_family','family_ratio'] else 'within','log2_mean':vals.mean(),'log2_q025':vals.quantile(.025),'log2_q975':vals.quantile(.975),'geometric_mean':2**vals.mean(),'draws':1000})
q=pd.DataFrame(qq);save(q,'figure2_components48.tsv');oldq=read(M/'tables/quantity_contrasts.tsv');oldq=oldq[oldq.scenario=='augmented_same_reference']
for r in q[q.version=='ordinary_normalized'].itertuples():
 oc='relative_family_abundance' if r.component=='family' else 'relative_loci_per_family';rr=oldq[(oldq.clade==r.clade)&(oldq.component==oc)].iloc[0];check('quantity_'+r.clade+r.component,np.isclose(r.geometric_mean,rr.geometric_mean))
save(pd.DataFrame(checks),'synchronization_checks.tsv');save(pd.DataFrame(manifest),'input_manifest.tsv')
(B/'SYNC_STATUS.json').write_text(json.dumps({'complete':True,'main_species':48,'matched_species_OGs':1940,'strict_species_OGs':1973,'strict_loci':2299,'original47_reproduced':True,'no_new_hypothesis_test':True,'new_resampling_seed':20260923,'TPM_unchanged':True},indent=2))
print(pd.DataFrame(alld).to_string(index=False));print(primary.to_string(index=False))
