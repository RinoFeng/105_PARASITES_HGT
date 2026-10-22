from pathlib import Path
import json
import numpy as np
import pandas as pd
O=Path(__file__).resolve().parents[1];CATS=list('EFGH');PRIMARY=['E','G','H'];SEED=20261003;R=5000
frames={label:pd.read_csv(O/'data'/fn,sep='\t') for label,fn in [('original48','original48_family_effects.tsv.gz'),('same_scope_old_background','family_effects_old_fixed_scope.tsv.gz'),('same_scope_controlled_background','family_effects_controlled_fixed_scope.tsv.gz')]}
OGS=sorted(set().union(*(set(x.OG) for x in frames.values())));oi={o:i for i,o in enumerate(OGS)};curves=[];metrics_rows=[];mass_dict={};positive_cov={}
def reweight(x):
 x=x.copy();nf=x.groupby('species').size();ns=x.groupby('clade').species.nunique();x['weight']=1/x.species.map(nf)/x.clade.map(ns)/len(ns);return x

def metric(v):
 z=np.sort(np.asarray(v)[np.asarray(v)>1e-12])[::-1]
 if not len(z):return {'positive_OGs':0,'top5_share':np.nan,'n50':np.nan,'effective_OGs':np.nan,'normalized_HHI':np.nan,'gross_positive':0}
 w=z/z.sum();h=(w*w).sum();return {'positive_OGs':len(z),'top5_share':w[:5].sum(),'n50':int(np.searchsorted(np.cumsum(w),.5)+1),'effective_OGs':1/h,'normalized_HHI':(h-1/len(z))/(1-1/len(z)) if len(z)>1 else np.nan,'gross_positive':z.sum()}

def masses(x):
 out={}
 ix=x.OG.map(oi).to_numpy()
 for c in CATS:out[c]=np.bincount(ix,weights=x.weight*np.maximum(x[c],0),minlength=len(OGS))
 return out

def contrast(d):return d['H']-.5*(d['E']+d['G'])
for label,x in frames.items():
 mass_dict[label]=masses(x)
 for cat in CATS:
  m=mass_dict[label][cat];met=metric(m);metrics_rows.append({'scenario':label,'category':cat,**met,'gross_negative':float((x.weight*np.minimum(x[cat],0)).sum()),'net_effect':float((x.weight*x[cat]).sum()),'species':x.species.nunique(),'species_OGs':len(x)})
  pos=x[x[cat]>1e-12]
  cov=pos.groupby('OG').agg(positive_species=('species','nunique'),positive_clades=('clade','nunique'),positive_loci=('loci','sum')).reset_index();cov['category']=cat;cov['scenario']=label;cov['positive_contribution']=cov.OG.map({g:m[oi[g]] for g in cov.OG});positive_cov[(label,cat)]=cov
  for idx in np.flatnonzero(m>1e-12):curves.append({'scenario':label,'category':cat,'OG':OGS[idx],'positive_contribution':m[idx]})
met=pd.DataFrame(metrics_rows);met.to_csv(O/'results/concentration_scenarios.tsv',sep='\t',index=False);pd.DataFrame(curves).to_csv(O/'results/OG_positive_contributions.tsv',sep='\t',index=False);pd.concat(positive_cov.values()).to_csv(O/'results/OG_distribution_opportunity.tsv',sep='\t',index=False)
expected={'E':(79,11,.329703),'F':(42,7,.411574),'G':(67,14,.216054),'H':(46,5,.514339)}
for c,(n,n50,t5) in expected.items():
 z=met[(met.scenario=='original48')&(met.category==c)].iloc[0];assert z.positive_OGs==n and z.n50==n50 and abs(z.top5_share-t5)<1e-6
# Fixed frozen original top-five lists, no category-dependent taxon deletion.
a0=frames['original48'];defs={c:[OGS[i] for i in np.argsort(-mass_dict['original48'][c])[:5]] for c in CATS};pd.DataFrame([{'category':c,'OGs':';'.join(gs)} for c,gs in defs.items()]).to_csv(O/'results/frozen_top5_OGs.tsv',sep='\t',index=False)
delrows=[]
for label,x in frames.items():
 for c in CATS:
  rem=x[~x.OG.isin(defs[c])].copy();valid=set(rem.species)==set(x.species)
  if valid:rem=reweight(rem);base=(x.weight*x[c]).sum();after=(rem.weight*rem[c]).sum()
  else:base=after=np.nan
  delrows.append({'scenario':label,'category':c,'species_before':x.species.nunique(),'species_after':rem.species.nunique(),'minimum_families_after':int(rem.groupby('species').size().min()),'same_species':valid,'baseline_pp':100*base,'remaining_pp':100*after,'proportional_reduction':(base-after)/base if base>0 else np.nan})
pd.DataFrame(delrows).to_csv(O/'results/deletion_fixed_species.tsv',sep='\t',index=False)
# Fixed scope, leave one clade; this is sensitivity, not new independent replication.
clrows=[]
x=frames['same_scope_controlled_background']
for scope,xx in [('all_four',x)]+[('only_'+cl,x[x.clade==cl]) for cl in sorted(x.clade.unique())]+[('omit_'+cl,x[x.clade!=cl]) for cl in sorted(x.clade.unique())]:
 xx=reweight(xx);mm=masses(xx);zz={c:metric(mm[c]) for c in CATS};th=contrast({c:zz[c]['top5_share'] for c in CATS})
 for c in CATS:clrows.append({'scope':scope,'category':c,'species':xx.species.nunique(),**zz[c],'H_minus_mean_EG':th})
pd.DataFrame(clrows).to_csv(O/'results/clade_sensitivity.tsv',sep='\t',index=False)
# Identical OG priorities preserve cross-category identity in finite subset draws.
priority=np.random.default_rng(SEED).random((R,len(OGS)));draws=[]
primary_label='same_scope_controlled_background';primary_m=int(min(metric(mass_dict[primary_label][c])['positive_OGs'] for c in PRIMARY));sizes=sorted(set([20,30,primary_m]));rare_stats=[]
def sample_stats(m,selected):
 vals=m[selected];tot=vals.sum(axis=1);sort=np.sort(vals,axis=1)[:,::-1];w=sort/tot[:,None];n50=(np.cumsum(w,axis=1)<.5).sum(axis=1)+1
 return {'top5_share':w[:,:5].sum(axis=1),'n50':n50,'effective_OGs':1/(w*w).sum(axis=1),'retained_mass_fraction':tot/m.sum()}
for label,mm in mass_dict.items():
 for size in sizes:
  ans={}
  for c in CATS:
   ids=np.flatnonzero(mm[c]>1e-12)
   if len(ids)<size:continue
   selected=ids[np.argpartition(priority[:,ids],size-1,axis=1)[:,:size]];stats=sample_stats(mm[c],selected);ans[c]=stats['top5_share']
   for stat,v in stats.items():rare_stats.append({'scenario':label,'design':'equal_positive_OG_count','m':size,'category':c,'statistic':stat,'mean':np.mean(v),'q025':np.quantile(v,.025),'q975':np.quantile(v,.975),'draws':R})
  if all(c in ans for c in PRIMARY):
   theta=contrast(ans)
   for j in range(R):draws.append({'scenario':label,'design':'equal_positive_OG_count','m':size,'replicate':j+1,'H_minus_mean_EG':theta[j],**{c:ans[c][j] if c in ans else np.nan for c in CATS}})
# Shared distribution strata are determined without using contribution magnitudes.
layer={};scopes=[];common_bins={}
for c in PRIMARY:
 z=positive_cov[(primary_label,c)].copy();z['species_bin']=np.where(z.positive_species==1,'1',np.where(z.positive_species<=3,'2-3','4plus'));z['clade_bin']=np.where(z.positive_clades==1,'1','2plus');z['stratum']=z.species_bin+'__'+z.clade_bin;layer[c]=z
for key in sorted(set().union(*(set(z.stratum) for z in layer.values()))):
 counts={c:int((layer[c].stratum==key).sum()) for c in PRIMARY};common_bins[key]=min(counts.values());scopes.append({'stratum':key,**counts,'sample_each_category':common_bins[key]})
pd.DataFrame(scopes).to_csv(O/'results/common_distribution_strata.tsv',sep='\t',index=False)
ms=sum(common_bins.values());ans={}
for c in PRIMARY:
 z=layer[c];selected=[];eligible=[]
 for key,n in common_bins.items():
  if n==0:continue
  ids=z.loc[z.stratum==key,'OG'].map(oi).to_numpy();eligible.extend(ids.tolist());selected.append(ids[np.argpartition(priority[:,ids],n-1,axis=1)[:,:n]])
 if ms<5:continue
 selected=np.concatenate(selected,axis=1);m=mass_dict[primary_label][c];stats=sample_stats(m,selected);ans[c]=stats['top5_share']
 for stat,v in stats.items():rare_stats.append({'scenario':primary_label,'design':'equal_count_and_distribution','m':ms,'category':c,'statistic':stat,'mean':np.mean(v),'q025':np.quantile(v,.025),'q975':np.quantile(v,.975),'draws':R,'eligible_pool_mass_fraction':m[np.unique(eligible)].sum()/m.sum()})
if all(c in ans for c in PRIMARY):
 theta=contrast(ans)
 for j in range(R):draws.append({'scenario':primary_label,'design':'equal_count_and_distribution','m':ms,'replicate':j+1,'H_minus_mean_EG':theta[j],**{c:ans[c][j] for c in PRIMARY}})
draw=pd.DataFrame(draws);draw.to_csv(O/'data/finite_subset_draws.tsv.gz',sep='\t',index=False);pd.DataFrame(rare_stats).to_csv(O/'results/finite_subset_metric_summary.tsv',sep='\t',index=False)
summary=[]
for (sc,de,sz),z in draw.groupby(['scenario','design','m']):
 v=z.H_minus_mean_EG.to_numpy();summary.append({'scenario':sc,'design':de,'m':sz,'mean_difference_pp':100*np.mean(v),'q025_pp':100*np.quantile(v,.025),'q975_pp':100*np.quantile(v,.975),'fraction_above_zero':np.mean(v>0),'draws':len(v),'interval_type':'finite_subset_variation_not_sampling_CI'})
pd.DataFrame(summary).to_csv(O/'results/finite_subset_contrast_summary.tsv',sep='\t',index=False)
# Matching quality, independently of functional concentration.
bal=pd.read_csv(O/'data/covariate_pair_differences.tsv.gz',sep='\t');br=[]
for label,z in bal.groupby('scenario'):
 for col in ['length_log_ratio','member_log_ratio']+[c+'_difference' for c in ['eggnog_fraction','informative_cog_fraction','ko_fraction','pfam_fraction']]:
  br.append({'scenario':label,'covariate':col,'weighted_mean_difference':np.average(z[col],weights=z.pair_weight),'weighted_mean_absolute_difference':np.average(abs(z[col]),weights=z.pair_weight)})
 br.append({'scenario':label,'covariate':'member_bin_mismatch','weighted_mean_difference':np.average(~z.same_member_bin,weights=z.pair_weight),'weighted_mean_absolute_difference':np.average(~z.same_member_bin,weights=z.pair_weight)})
pd.DataFrame(br).to_csv(O/'results/matching_balance.tsv',sep='\t',index=False)
# Conditional ordinary-vs-ordinary construction: shared ordinary OG priorities across species/sets.
pairs=pd.read_csv(O/'data/matched_pairs_controlled_fixed_scope.tsv.gz',sep='\t');b=pd.read_csv(O/'data/controlled_background_profiles.tsv.gz',sep='\t').set_index(['species','OG']);bogs=sorted(pairs.background_OG.unique());bi={g:i for i,g in enumerate(bogs)};pdict={key:gg for key,gg in pairs.groupby(['species','candidate_OG'])};x=x.reset_index(drop=True);n=len(x)
ids=np.zeros((n,5),dtype=int);vals=np.zeros((n,5,4));mask=np.zeros((n,5),dtype=bool);counts=np.zeros(n,dtype=int)
for i,row in enumerate(x.itertuples()):
 ps=pdict[(row.species,row.OG)];counts[i]=len(ps)
 for j,og in enumerate(ps.background_OG):ids[i,j]=bi[og];vals[i,j]=b.loc[(row.species,og),CATS].to_numpy();mask[i,j]=True
assert counts.min()>=3
sums=vals.sum(axis=1);weight=x.weight.to_numpy();rng=np.random.default_rng(SEED+1);bgrows=[]
for rep in range(2000):
 pr=rng.random(len(bogs));chosen=np.argmin(np.where(mask,pr[ids],np.inf),axis=1);focus=vals[np.arange(n),chosen];contrast_vals=focus-(sums-focus)/(counts-1)[:,None];focus_ids=ids[np.arange(n),chosen];vp=rng.random(len(bogs));out={'replicate':rep+1};t5={}
 for ki,c in enumerate(CATS):
  bm=np.bincount(focus_ids,weights=weight*np.maximum(contrast_vals[:,ki],0),minlength=len(bogs));positive=np.flatnonzero(bm>1e-12);out[c+'_positive_OGs']=len(positive)
  if len(positive)>=20:
   take=positive[np.argpartition(vp[positive],19)[:20]];w=bm[take];t5[c]=np.sort(w)[-5:].sum()/w.sum();out[c+'_top5_at_m20']=t5[c]
  else:out[c+'_top5_at_m20']=np.nan
 out['main_valid']=all(c in t5 for c in PRIMARY);out['H_minus_mean_EG']=contrast(t5) if out['main_valid'] else np.nan;bgrows.append(out)
bgr=pd.DataFrame(bgrows);bgr.to_csv(O/'data/ordinary_reference_draws.tsv.gz',sep='\t',index=False);valid=bgr[bgr.main_valid];ob=draw[(draw.scenario==primary_label)&(draw.design=='equal_positive_OG_count')&(draw.m==20)];obs=float(ob.H_minus_mean_EG.mean());ref=valid.H_minus_mean_EG.to_numpy()
refsum={'reference_model':'ordinary_focus_vs_remaining_matched_ordinary;shared_OG_priority','draws':2000,'m':20,'valid_draws':len(valid),'valid_fraction':len(valid)/2000,'observed_mean_contrast_pp':100*obs,'observed_subset_q025_pp':100*ob.H_minus_mean_EG.quantile(.025),'observed_subset_q975_pp':100*ob.H_minus_mean_EG.quantile(.975),'reference_mean_pp':100*np.mean(ref) if len(ref) else None,'reference_q025_pp':100*np.quantile(ref,.025) if len(ref) else None,'reference_q975_pp':100*np.quantile(ref,.975) if len(ref) else None,'fraction_reference_ge_observed_mean':float(np.mean(ref>=obs)) if len(ref) else None,'formal_P_value':None,'interpretation':'conditional_background_diagnostic_not_HGT_label_randomization'}
refsum.update({'reference_adequacy_pass':len(valid)/2000>=.9,'tail_inference_supported':False,'decision':'insufficient_valid_fraction; do_not_interpret_tail_fraction_as_evidence_or_P' if len(valid)/2000<.9 else 'conditional_diagnostic_only; no_formal_P'})
(O/'results/ordinary_reference_summary.json').write_text(json.dumps(refsum,indent=2));(O/'results/analysis_summary.json').write_text(json.dumps({'status':'computed','seed':SEED,'finite_subset_draws':R,'primary_family_count':primary_m,'distribution_balanced_family_count':ms,'reference':refsum},indent=2));print(met[['scenario','category','positive_OGs','top5_share','n50']].to_string(index=False));print(pd.DataFrame(summary).to_string(index=False));print(json.dumps(refsum,indent=2))
