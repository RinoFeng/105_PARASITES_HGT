from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import solve_triangular
from audit import OUT, limited_tree, covariance

SEED=20260922
REPS=10000
ENDS={'Y':'Y_family_log2','D':'D_family_minus_within_log2','Z':'Z_loci_per_family_log2','T':'T_total_loci_log2'}

def setup(x,v):
    l=np.linalg.cholesky(v)
    w=solve_triangular(l,np.eye(len(v)),lower=True)
    wx=w@x
    inv=np.linalg.inv(wx.T@wx)
    a=inv@wx.T@w
    return w,wx,inv,a

def fit(x,y,v):
    w,wx,inv,a=setup(x,v)
    beta=a@y
    residual=w@(y-x@beta)
    df=len(y)-x.shape[1]
    s2=residual@residual/df
    se=np.sqrt(s2*inv[1,1])
    q=stats.t.ppf(.975,df); qsim=stats.t.ppf(1-.05/12,df)
    return dict(estimate=beta[1],se=se,df=df,t=beta[1]/se,p_raw=2*stats.t.sf(abs(beta[1]/se),df),ci95_low=beta[1]-q*se,ci95_high=beta[1]+q*se,ci_family_low=beta[1]-qsim*se,ci_family_high=beta[1]+qsim*se,residual_variance=s2), residual

def holm(p):
    p=np.asarray(p); order=np.argsort(p); a=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
    out=np.empty(len(p));out[order]=a
    return out

def wilson(k,n):
    p=k/n; z=stats.norm.ppf(.975); d=1+z*z/n
    mid=(p+z*z/(2*n))/d
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return mid-half,mid+half

def simulate(x,vfit,vtrue,rng):
    w,wx,inv,a=setup(x,vfit)
    y=np.linalg.cholesky(vtrue)@rng.standard_normal((len(x),REPS))
    beta=a@y
    residual=w@(y-x@beta)
    df=len(x)-x.shape[1]
    s2=(residual**2).sum(axis=0)/df
    t=beta[1]/np.sqrt(s2*inv[1,1])
    k=int(np.sum(np.abs(t)>stats.t.ppf(.975,df)))
    lo,hi=wilson(k,REPS)
    return dict(replicates=REPS,rejections=k,type_I=k/REPS,wilson95_low=lo,wilson95_high=hi,within_prespecified_calibration_range=.04<=k/REPS<=.06)

def recompute(frame,f,removed):
    ff=f[f.species.isin(frame.species)&~f.OG.isin(removed)]
    agg=ff.groupby(['species','role']).agg(families=('OG','size'),loci=('units','sum'))
    z=frame.copy()
    for col,role,measure in [('C','candidate','families'),('B','background','families'),('N','candidate','loci'),('M','background','loci')]:
        z[col]=[agg.loc[(sp,role),measure] for sp in z.species]
    z[ENDS['Y']]=np.log2(z.C/z.B)
    z[ENDS['Z']]=np.log2((z.N/z.C)/(z.M/z.B))
    z[ENDS['D']]=z[ENDS['Y']]-z[ENDS['Z']]
    z[ENDS['T']]=z[ENDS['Y']]+z[ENDS['Z']]
    return z

def main():
    protocol=OUT/'STATISTICAL_PROTOCOL.json'
    digest=hashlib.sha256(protocol.read_bytes()).hexdigest()
    s=pd.read_csv(OUT/'species_observation_audit.tsv',sep='\t')
    f=pd.read_csv(OUT/'family_inputs48.tsv.gz',sep='\t')
    results=[];calibration=[];deleted=[];shared=[];diagnostics=[];weightrows=[]
    for ci,clade in enumerate(['Apicomplexa','Euglenozoa']):
        z=s[s.clade.isin([clade,'Metazoa'])].copy().reset_index(drop=True)
        x=np.column_stack([np.ones(len(z)),z.clade.eq(clade).astype(float)])
        t=limited_tree(z)
        v=covariance(t,list(z.species_id))
        vp=covariance(limited_tree(z,'PMSF'),list(z.species_id))
        vg=covariance(t,list(z.species_id),scale='grafen')
        vl=covariance(t,list(z.species_id),lam=.5)
        variants={'LG_primary':v,'PMSF':vp,'Grafen':vg,'lambda0.5':vl,
                  'root0.25':covariance(limited_tree(z,fraction=.25),list(z.species_id)),
                  'root0.75':covariance(limited_tree(z,fraction=.75),list(z.species_id))}
        for name,vv in variants.items():
            for endpoint,col in ENDS.items():
                row,resid=fit(x,z[col].to_numpy(),vv)
                if endpoint in ['Z','T']:row['p_raw']=np.nan
                results.append(dict(comparison=clade,model=name,endpoint=endpoint,species=len(z),**row))
        # GLS weights: a clade contrast is not simply an equal-species mean.
        _,_,_,a=setup(x,v)
        for j,r in z.iterrows():weightrows.append(dict(comparison=clade,species=r.species,genus=r.genus,clade=r.clade,contrast_weight=a[1,j]))
        _,wx,_,_=setup(x,v)
        q=np.linalg.qr(wx,mode='complete')[0][:,2:]
        w=setup(x,v)[0]
        for endpoint in ['Y','D']:
            independent=q.T@w@z[ENDS[endpoint]].to_numpy()
            shap=stats.shapiro(independent)
            diagnostics.append(dict(comparison=clade,endpoint=endpoint,diagnostic='Shapiro_QR_residual_coordinates',statistic=shap.statistic,p_diagnostic_only=shap.pvalue,n=len(independent),note='coordinate-dependent diagnostic; not used to choose covariance or rescue P'))
        rng=np.random.default_rng(SEED+ci)
        # Matched fits verify implementation for every declared covariance.
        for name,vv in [('LG',v),('PMSF',vp),('Grafen',vg),('lambda0.5',vl)]:
            calibration.append(dict(comparison=clade,fitted_covariance=name,generating_covariance=name,kind='matched_model',**simulate(x,vv,vv,rng)))
        for name,vv in [('PMSF',vp),('Grafen',vg),('lambda0.5',vl)]:
            calibration.append(dict(comparison=clade,fitted_covariance='LG',generating_covariance=name,kind='misspecification_stress',**simulate(x,v,vv,rng)))
        for sd in [.5,2.]:
            mult=np.where(z.clade.eq('Metazoa'),sd,1.)
            vv=v*mult[:,None]*mult[None,:]
            calibration.append(dict(comparison=clade,fitted_covariance='LG',generating_covariance=f'nematode_SD_x{sd}',kind='heteroscedastic_stress',**simulate(x,v,vv,rng)))
        for genus in z.genus.unique():
            keep=z.genus.ne(genus).to_numpy()
            for endpoint in ['Y','D']:
                row,_=fit(x[keep],z.loc[keep,ENDS[endpoint]].to_numpy(),v[np.ix_(keep,keep)])
                deleted.append(dict(comparison=clade,omitted_genus=genus,endpoint=endpoint,species=int(keep.sum()),estimate=row['estimate'],ci95_low=row['ci95_low'],ci95_high=row['ci95_high']))
        keep=z.species.ne('Angiostrongylus costaricensis').to_numpy()
        for endpoint,col in ENDS.items():
            row,_=fit(x[keep],z.loc[keep,col].to_numpy(),v[np.ix_(keep,keep)])
            if endpoint in ['Z','T']:row['p_raw']=np.nan
            results.append(dict(comparison=clade,model='original47_without_Angiostrongylus',endpoint=endpoint,species=int(keep.sum()),**row))
        occupancy=f[(f.species.isin(z.species))&(f.role=='candidate')].groupby('OG').species.nunique().reset_index(name='species_occupancy')
        top=occupancy.sort_values(['species_occupancy','OG'],ascending=[False,True]).head(5)
        for _,r in top.iterrows():shared.append(dict(comparison=clade,OG=r.OG,species_occupancy=r.species_occupancy))
        zz=recompute(z,f,set(top.OG))
        zz.to_csv(OUT/f'{clade}_after_sharedOG_deletion.tsv',sep='\t',index=False)
        for endpoint,col in ENDS.items():
            row,_=fit(x,zz[col].to_numpy(),v)
            if endpoint in ['Z','T']:row['p_raw']=np.nan
            results.append(dict(comparison=clade,model='remove_top5_shared_candidate_OG',endpoint=endpoint,species=len(z),**row))
    out=pd.DataFrame(results)
    for name,ix in out[out.endpoint.isin(['Y','D'])].groupby('model').groups.items():
        vals=out.loc[ix,'p_raw'].to_list()+[1,1]
        assert len(vals)==6
        out.loc[ix,'p_Holm6']=holm(vals)[:4]
    primary=out[(out.model=='LG_primary')&out.endpoint.isin(['Y','D'])].copy()
    primary['effect_ratio']=2**primary.estimate
    primary['ratio_ci95_low']=2**primary.ci95_low;primary['ratio_ci95_high']=2**primary.ci95_high
    primary['inference_status']='conditional_working_covariance_only_not_ascertainment_adjusted'
    primary.to_csv(OUT/'primary_conditional_tests.tsv',sep='\t',index=False)
    out.to_csv(OUT/'all_model_results.tsv',sep='\t',index=False)
    pd.DataFrame(calibration).to_csv(OUT/'null_calibration.tsv',sep='\t',index=False)
    pd.DataFrame(deleted).to_csv(OUT/'leave_one_genus_effects.tsv',sep='\t',index=False)
    pd.DataFrame(shared).to_csv(OUT/'shared_OG_deletion_definition.tsv',sep='\t',index=False)
    pd.DataFrame(diagnostics).to_csv(OUT/'residual_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(weightrows).to_csv(OUT/'GLS_contrast_weights.tsv',sep='\t',index=False)
    # Same descriptive genus scale for all three lineages; no hypothesis tests.
    gen=s.groupby(['clade','genus'])[[*ENDS.values()]].mean().reset_index()
    desc=[]
    for clade in ['Apicomplexa','Euglenozoa','Amoebozoa']:
        for ep,col in ENDS.items():
            xx=gen[gen.clade==clade][col]; yy=gen[gen.clade=='Metazoa'][col]
            desc.append(dict(comparison=clade,endpoint=ep,focal_genera=len(xx),nematode_genera=len(yy),equal_genus_log2_difference=xx.mean()-yy.mean(),equal_genus_ratio=2**(xx.mean()-yy.mean()),minimum_pair_difference=min(xx)-max(yy),maximum_pair_difference=max(xx)-min(yy)))
    pd.DataFrame(desc).to_csv(OUT/'equal_genus_description.tsv',sep='\t',index=False)
    checks=[]
    for (clade,model),g in out.groupby(['comparison','model']):
        est=g.set_index('endpoint').estimate
        checks.append(dict(check=f'{clade}/{model}/T_equals_Y_plus_Z',passed=bool(np.isclose(est['T'],est['Y']+est['Z']))))
        checks.append(dict(check=f'{clade}/{model}/D_equals_Y_minus_Z',passed=bool(np.isclose(est['D'],est['Y']-est['Z']))))
    for clade in ['Apicomplexa','Euglenozoa']:
        g=out[out.comparison==clade].set_index(['model','endpoint'])
        for model in ['root0.25','root0.75']:
            for ep in ENDS:
                checks.append(dict(check=f'{clade}/{model}/{ep}/root_invariance',passed=bool(np.allclose(g.loc[('LG_primary',ep),['estimate','se']].to_numpy(float),g.loc[(model,ep),['estimate','se']].to_numpy(float)))))
    checks.append(dict(check='protocol_unchanged',passed=hashlib.sha256(protocol.read_bytes()).hexdigest()==digest))
    pd.DataFrame(checks).to_csv(OUT/'verification_checks.tsv',sep='\t',index=False)
    assert all(x['passed'] for x in checks)
    (OUT/'ANALYSIS_STATUS.json').write_text(json.dumps({'completed':True,'protocol_sha256':digest,'seed':SEED,'simulations_per_scenario':REPS,'Amoebozoa_tests':'not_run_single_genus','logs_read':False,'safety_skipped_tasks':[]},indent=2))
    print(primary.to_string(index=False))
    print(pd.DataFrame(calibration).to_string(index=False))
    print(pd.DataFrame(diagnostics).to_string(index=False))

if __name__=='__main__': main()
