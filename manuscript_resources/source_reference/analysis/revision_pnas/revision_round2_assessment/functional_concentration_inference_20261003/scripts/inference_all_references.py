from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.special import softmax
O=Path(__file__).resolve().parents[1];SEED=202610034;R=9999;CAL=1000;CATS=list('EFGH');PRI=[0,2,3]
x=pd.read_csv(O/'data/common_analysis_profiles.tsv.gz',sep='\t');roles=['candidate','reference','diagnostic_A','diagnostic_B'];blocks=sorted(x.block.unique());rng=np.random.default_rng(SEED)
M=[]
for role in roles:
 z=x[x.analysis_role==role].copy();z[CATS]=z[CATS].multiply(z.weight,axis=0);M.append(z.groupby('block')[CATS].sum().reindex(blocks).to_numpy())
cube=np.stack(M,axis=1);n=len(blocks)
def calc(data,choices):
 N,K,_=data.shape;r=len(choices);A=data[np.arange(N)[None,:],choices,:];B=np.broadcast_to(data,(r,)+data.shape).copy();B[np.arange(r)[:,None],np.arange(N)[None,:],choices,:]=0;B=B.reshape(r,N*K,4)
 sa=A.sum(axis=1);sb=B.sum(axis=1)
 with np.errstate(invalid='ignore',divide='ignore'):
  ta=np.partition(A,max(0,N-5),axis=1)[:,-5:,:].sum(axis=1)/sa;tb=np.partition(B,max(0,N*K-5),axis=1)[:,-5:,:].sum(axis=1)/sb
 d=ta-tb;D=d[:,3]-.5*(d[:,0]+d[:,2]);return D,ta,tb

def calc_chunks(data,choices):
 return np.concatenate([calc(data,choices[i:i+256])[0] for i in range(0,len(choices),256)])

def dp_build(signatures,target,prob):
 N,K,_=signatures.shape;shape=tuple((target+1).tolist());dp=np.zeros((N+1,)+shape);dp[N,0,0,0]=1
 for i in range(N-1,-1,-1):
  for k in range(K):
   v=signatures[i,k];dest=tuple(slice(int(d),int(t)+1) for d,t in zip(v,target));src=tuple(slice(0,int(t-d)+1) for d,t in zip(v,target))
   if np.any(v>target):continue
   dp[(i,)+dest]+=prob[i,k]*dp[(i+1,)+src]
 return dp

def sample(data,N,cond,rng,prob=None):
 g,k,_=data.shape
 if prob is None:prob=np.full((g,k),1/k)
 sig=(data[:,:,PRI]>1e-14).astype(np.int16);target=sig[:,0,:].sum(axis=0);total=sig.sum(axis=(0,1))
 if cond:
  dp=dp_build(sig,target,prob);mass=float(dp[(0,)+tuple(target)])
  assert mass>0,'zero conditional probability'
  choices=np.zeros((N,g),dtype=np.int8);remain=np.tile(target,(N,1))
  for i in range(g):
   probs=np.zeros((N,k))
   for j in range(k):
    rr=remain-sig[i,j];ok=np.all(rr>=0,axis=1)
    if ok.any():zz=rr[ok];probs[ok,j]=prob[i,j]*dp[i+1,zz[:,0],zz[:,1],zz[:,2]]
   den=probs.sum(axis=1);assert np.all(den>0);probs/=den[:,None];sel=(rng.random(N)[:,None]>np.cumsum(probs,axis=1)).sum(axis=1);sel=np.minimum(sel,k-1);choices[:,i]=sel;remain-=sig[i,sel]
  assert np.all(remain==0)
  return choices,{'draws':N,'conditional_event_probability':mass,'sampling':'exact_dynamic_program_conditional_sampling','target_counts':target.tolist()}
 else:
  pieces=[];got=0;proposed=0;accepted=0
  while got<N:
   size=max(1024,N-got);u=rng.random((size,g));sel=(u[:,:,None]>np.cumsum(prob,axis=1)[None,:,:]).sum(axis=2);sel=np.minimum(sel,k-1).astype(np.int8)
   cnt=sig[np.arange(g)[None,:],sel,:].sum(axis=1);ok=np.all((cnt>0)&(cnt<total[None,:]),axis=1);proposed+=size;accepted+=int(ok.sum());take=sel[ok][:N-got];pieces.append(take);got+=len(take)
   if proposed>20000000:raise RuntimeError('too few admissible assignments')
  return np.concatenate(pieces),{'draws':N,'proposed':proposed,'admissible_fraction':accepted/proposed,'sampling':'uniform_assignments_conditional_on_nonzero_denominators','target_counts':None}

def twop(obs,ref):
 s=np.sort(ref);v=np.asarray(obs);N=len(s);left=(np.searchsorted(s,v+1e-12,side='right')+1)/(N+1);right=(N-np.searchsorted(s,v-1e-12,side='left')+1)/(N+1);return np.minimum(1,2*np.minimum(left,right))
def wilson(k,n):
 z=1.95996398454;p=k/n;den=1+z*z/n;ctr=(p+z*z/(2*n))/den;rad=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den;return ctr-rad,ctr+rad
# Dynamic program probabilities checked by full enumeration on a small synthetic quartet example.
import itertools
sigtest=np.array([[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],[[1,0,0],[0,0,0],[1,1,0],[0,0,1]],[[0,1,0],[0,0,0],[1,0,1],[0,0,1]],[[0,0,1],[0,1,0],[0,0,0],[1,0,0]]],dtype=np.int16)
targ=np.array([1,1,1]);pt=softmax(np.arange(16).reshape(4,4)/20,axis=1);dp=dp_build(sigtest,targ,pt);exact=0
for labels in itertools.product(range(4),repeat=4):
 if np.array_equal(sigtest[np.arange(4),labels].sum(axis=0),targ):exact+=np.prod(pt[np.arange(4),labels])
assert abs(dp[0,1,1,1]-exact)<1e-12
observed,ta,tb=calc(cube,np.zeros((1,n),dtype=np.int8));observed=float(observed[0]);assert np.isfinite(observed)
rows=[]
for role,m,vals in [('candidate',cube[:,0,:],ta[0]),('all_three_ordinary',cube[:,1:,:].reshape(n*3,4),tb[0])]:
 for j,c in enumerate(CATS):rows.append({'role':role,'category':c,'positive_OGs':int((m[:,j]>1e-14).sum()),'top5_share':float(vals[j]),'function_mass':float(m[:,j].sum()/(3 if role=='all_three_ordinary' else 1))})
pd.DataFrame(rows).to_csv(O/'results/all_references_observed.tsv',sep='\t',index=False)
# Pooled background covariate balance; quartet residual assignment probabilities for sensitivity.
ca=x[x.analysis_role=='candidate'].set_index(['block','species']);balance=[]
for name in ['length','units','eggnog_fraction','informative_cog_fraction','ko_fraction','pfam_fraction']:
 aa=ca[name].to_numpy();w=ca.weight.to_numpy();zz=x[x.analysis_role!='candidate'];bb=zz[name].to_numpy();wb=zz.weight.to_numpy()/3
 if name in ['length','units']:aa=np.log(aa);bb=np.log(bb)
 ma=np.average(aa,weights=w);mb=np.average(bb,weights=wb);sd=np.sqrt((np.average((aa-ma)**2,weights=w)+np.average((bb-mb)**2,weights=wb))/2);smd=(ma-mb)/sd if sd else 0
 balance.append({'covariate':name,'scale':'log' if name in ['length','units'] else 'fraction','candidate_mean':ma,'ordinary_mean':mb,'SMD':smd})
pd.DataFrame(balance).to_csv(O/'results/all_references_balance.tsv',sep='\t',index=False)
score=np.zeros((n,4))
for i,b in enumerate(blocks):
 for j,role in enumerate(roles):
  z=x[(x.block==b)&(x.analysis_role==role)];score[i,j]=.5*(np.average(np.log(z.length),weights=z.weight)/np.log(1.25)+np.average(np.log(z.units),weights=z.weight)/np.log(2))
 score[i]-=(score[i].min()+score[i].max())/2;score[i]/=max(1,np.ptp(score[i]))
results=[];cal=[];sens=[];drawrows=[];sampling=[];ordinary=[]
for cond in [False,True]:
 label='fixed_total_family_ratio_1to3' if not cond else 'also_fixed_EGH_family_counts'
 choices,info=sample(cube,R+CAL,cond,rng);values=calc_chunks(cube,choices);ref=values[:R];probe=values[R:];assert np.isfinite(values).all();p=float(twop(observed,ref));variable=np.ptp(ref)>1e-12
 info.update({'design':label,'unique_assignments':len(np.unique(choices,axis=0)),'distinct_statistics_rounded12':len(np.unique(np.round(ref,12)))});sampling.append(info)
 cps=twop(probe,ref);k=int((cps<=.05).sum());lo,hi=wilson(k,CAL);cal.append({'design':label,'replicates':CAL,'rejections_005':k,'rate':k/CAL,'wilson_lower':lo,'wilson_upper':hi,'inflation_detected':lo>.05})
 results.append({'design':label,'D_pp':100*observed,'null_mean_pp':float(100*np.mean(ref)),'null_q025_pp':float(100*np.quantile(ref,.025)),'null_q975_pp':float(100*np.quantile(ref,.975)),'P_two_sided_conditional':p,'statistic_variable':variable,'uniform_calibration_inflation':lo>.05,'status':'exploratory_conditional_test'})
 for i,v in enumerate(values):drawrows.append({'design':label,'replicate':i+1,'purpose':'reference' if i<R else 'calibration','D':v})
 print(label,'D_pp',100*observed,'P',p,'calibration',k/CAL,flush=True)
 for gamma in [1.25,1.5,2.0]:
  for direction in [-1,1]:
   pr=softmax(direction*np.log(gamma)*score,axis=1);ch,inf=sample(cube,CAL,cond,rng,prob=pr);vv=calc_chunks(cube,ch);pp=twop(vv,ref);k2=int((pp<=.05).sum());l2,h2=wilson(k2,CAL)
   sens.append({'design':label,'gamma_bound':gamma,'direction':direction,'replicates':CAL,'largest_within_group_probability_ratio':float((pr.max(axis=1)/pr.min(axis=1)).max()),'rejection_rate':k2/CAL,'wilson_lower':l2,'wilson_upper':h2,'inflation_detected':l2>.05})
# Ordinary-only 1:2 diagnostic: focus is the previously reserved diagnostic_A; no role chosen by effect.
bgcube=cube[:,[2,1,3],:];bo,_,_=calc(bgcube,np.zeros((1,n),dtype=np.int8));bo=float(bo[0])
for cond in [False,True]:
 label='fixed_total_family_ratio_1to2' if not cond else 'also_fixed_EGH_family_counts'
 if np.isfinite(bo):
  ch,inf=sample(bgcube,R,cond,rng);v=calc_chunks(bgcube,ch);ordinary.append({'design':label,'D_pp':100*bo,'P_diagnostic_only':float(twop(bo,v)),'null_q025_pp':100*np.quantile(v,.025),'null_q975_pp':100*np.quantile(v,.975),'not_same_reference_ratio_as_main':True})
 else:ordinary.append({'design':label,'status':'undefined','not_same_reference_ratio_as_main':True})
tests=pd.DataFrame(results);p=tests.P_two_sided_conditional.to_numpy();adj=np.empty(2);order=np.argsort(p);last=0
for rank,i in enumerate(order):last=max(last,min(1,(2-rank)*p[i]));adj[i]=last
tests['P_Holm_two_tests']=adj;tests.to_csv(O/'results/all_references_tests.tsv',sep='\t',index=False)
pd.DataFrame(cal).to_csv(O/'results/all_references_calibration.tsv',sep='\t',index=False);pd.DataFrame(sens).to_csv(O/'results/all_references_assignment_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(ordinary).to_csv(O/'results/all_references_ordinary_diagnostic.tsv',sep='\t',index=False);pd.DataFrame(sampling).to_csv(O/'results/all_references_sampling.tsv',sep='\t',index=False);pd.DataFrame(drawrows).to_csv(O/'data/all_references_randomization_draws.tsv.gz',sep='\t',index=False)
summary={'status':'supplementary_conditional_inference_computed','original_1to1_test_status':'undefined_E_reference_denominator','quartets':n,'candidate_OGs':n,'ordinary_OGs':3*n,'species':x.species.nunique(),'seed':SEED,'D_pp':100*observed,'tests':json.loads(tests.to_json(orient='records')),'maximum_absolute_SMD':float(np.max(np.abs(pd.DataFrame(balance).SMD))),'dp_probability_enumeration_check':True,'scope':'fixed_common_observed_subset;not_48_species_inference','interpretation':'conditional_on_quartet_label_exchangeability;exploratory_after_denominator_audit'}
(O/'results/all_references_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
