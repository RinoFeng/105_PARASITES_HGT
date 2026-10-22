from pathlib import Path
from collections import defaultdict
import numpy as np,pandas as pd,json,math
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
O=Path(__file__).resolve().parents[1];R=O.parent;SEED=202610032
f=pd.read_csv(R/'metazoa_extension_20260922/phylogenetic_comparison/family_inputs48.tsv.gz',sep='\t');a=pd.read_csv(R/'manuscript_refinement_20260923/tables/family_effects_augmented_same_reference.tsv.gz',sep='\t');target=sorted(a.OG.unique());covs=['eggnog_fraction','informative_cog_fraction','ko_fraction','pfam_fraction'];sps=sorted(f.species.unique());si={s:i for i,s in enumerate(sps)}
def abin(v):return 0 if v==0 else 3 if v==1 else 1 if v<=.5 else 2
def ubin(v):return 0 if v==1 else 1 if v==2 else 2 if v<=5 else 3
f=f.sort_values(['OG','species']);families={g:z.set_index('species') for g,z in f.groupby('OG')};bits={g:sum(1<<si[s] for s in z.index) for g,z in families.items()};bg=f[f.role=='background'].copy();bgs=sorted(bg.OG.unique());bi={b:i for i,b in enumerate(bgs)}
bg['key']=[(r.species,ubin(r.units))+tuple(abin(getattr(r,c)) for c in covs) for r in bg.itertuples()];pools={k:z.sort_values('length') for k,z in bg.groupby('key')};edges={};audit=[]
for idx,g in enumerate(target):
 x=families[g];compatible=defaultdict(list)
 for sp,r in x.iterrows():
  key=(sp,ubin(r.units))+tuple(abin(r[c]) for c in covs);z=pools.get(key)
  if z is None:continue
  lens=z.length.to_numpy();lo=np.searchsorted(lens,.8*r.length,side='left');hi=np.searchsorted(lens,1.25*r.length,side='right')
  for b in z.OG.iloc[lo:hi]:compatible[b].append(sp)
 selected={}
 for b,ss in compatible.items():
  if len(ss)<math.ceil(.8*len(x)):continue
  jac=bin(bits[g]&bits[b]).count('1')/bin(bits[g]|bits[b]).count('1')
  if jac<.5:continue
  cost=float(np.mean(abs(np.log(families[b].loc[ss,'length'].to_numpy()/x.loc[ss,'length'].to_numpy()))))
  selected[b]={'species':ss,'cost':cost,'jaccard':jac}
 edges[g]=selected;audit.append({'candidate_OG':g,'eligible_ordinary_OGs':len(selected),'original_species':len(x)})
# Sparse global non-reuse matching with a unique dummy for each slot.
rr=[];cc=[];vv=[];nr=3*len(target);nb=len(bgs)
for i,g in enumerate(target):
 for k in range(3):
  slot=i*3+k
  for b,e in edges[g].items():rr.append(slot);cc.append(bi[b]);vv.append(1+e['cost']+bi[b]*1e-12)
  rr.append(slot);cc.append(nb+slot);vv.append(10000.)
mat=coo_matrix((vv,(rr,cc)),shape=(nr,nb+nr)).tocsr();ridx,cidx=min_weight_full_bipartite_matching(mat);alloc=defaultdict(list)
for i,j in zip(ridx,cidx):
 if j<nb:alloc[target[i//3]].append(bgs[j])
records=[];rng=np.random.default_rng(SEED)
for ar in audit:
 g=ar['candidate_OG'];refs=alloc[g];ar['allocated']=len(refs);ar['status']='fewer_than_three_disjoint_references'
 if len(refs)!=3:continue
 ss=set.intersection(*(set(edges[g][b]['species']) for b in refs));ss=sorted(ss);ar['common_species']=len(ss)
 if len(ss)<math.ceil(.8*len(families[g])):ar['status']='common_coverage_below_80pct';continue
 lens=np.vstack([families[t].loc[ss,'length'].to_numpy() for t in [g]+refs])
 if np.max(lens.max(axis=0)/lens.min(axis=0))>1.25+1e-12:ar['status']='quartet_length_span';continue
 refs=list(np.array(sorted(refs))[rng.permutation(3)]);ar['status']='retained'
 records.append({'block':len(records)+1,'candidate_OG':g,'reference_OG':refs[0],'diagnostic_A_OG':refs[1],'diagnostic_B_OG':refs[2],'common_species':';'.join(ss),'candidate_original_species':';'.join(families[g].index),'reference_original_species':';'.join(families[refs[0]].index),'minimum_jaccard':min(edges[g][b]['jaccard'] for b in refs)})
q=pd.DataFrame(records);pd.DataFrame(audit).to_csv(O/'results/common_matching_audit.tsv',sep='\t',index=False);q.to_csv(O/'data/common_quartets_before_scope.tsv',sep='\t',index=False)
counts=defaultdict(int)
for ss in q.common_species:
 for sp in ss.split(';'):counts[sp]+=1
scope=f[['species','clade']].drop_duplicates().sort_values('species');scope['quartets']=scope.species.map(counts).fillna(0).astype(int);scope['retained']=scope.quartets>=10;scope.to_csv(O/'results/common_species_scope.tsv',sep='\t',index=False);keep=set(scope[scope.retained].species)
q=q[q.common_species.map(lambda x:bool(set(x.split(';'))&keep))].copy();q['analysis_species']=q.common_species.map(lambda x:';'.join(sorted(set(x.split(';'))&keep)));q.to_csv(O/'data/common_analysis_quartets.tsv',sep='\t',index=False)
rows=[]
for r in q.itertuples():
 for role,g in [('candidate',r.candidate_OG),('reference',r.reference_OG),('diagnostic_A',r.diagnostic_A_OG),('diagnostic_B',r.diagnostic_B_OG)]:
  z=families[g].loc[r.analysis_species.split(';')].reset_index().copy();z['analysis_role']=role;z['block']=r.block;rows.append(z)
out=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
if len(out):
 nf=out[out.analysis_role=='candidate'].groupby('species').size();ns=scope[scope.retained].groupby('clade').size();out['weight']=1/out.species.map(nf)/out.clade.map(ns)/len(ns);out.to_csv(O/'data/common_analysis_profiles.tsv.gz',sep='\t',index=False)
 assert all(np.isclose(w,1) for w in out.groupby('analysis_role').weight.sum())
 assert len(set(q.reference_OG)|set(q.diagnostic_A_OG)|set(q.diagnostic_B_OG))==3*len(q)
summary={'status':'common_design_ready' if len(q) else 'no_common_scope','target_OGs':len(target),'quartets_before_scope':len(records),'analysis_quartets':len(q),'species':len(keep),'clades':scope[scope.retained].groupby('clade').size().to_dict(),'species_candidate_OG_units':len(out)//4,'matching_reason_counts':pd.DataFrame(audit).status.value_counts().to_dict(),'seed':SEED,'function_values_used_for_selection':False}
(O/'results/common_design_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
