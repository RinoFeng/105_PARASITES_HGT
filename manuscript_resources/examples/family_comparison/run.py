#!/usr/bin/env python3
"""Reproduce fixed-panel quantity summaries from published species counts.
This is a downstream reproduction example, not an HGT discovery pipeline.
"""
from pathlib import Path
import argparse,json,platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
CLADES=['Amoebozoa','Apicomplexa','Euglenozoa','Metazoa']
COMPONENTS=['relative_family_abundance','relative_loci_per_family','relative_locus_abundance']

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=HERE/'output')
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(HERE/'input/species_counts.tsv',sep='\t')
    assert d.species.is_unique and len(d)==48
    assert d.groupby('clade').size().to_dict()==dict(Amoebozoa=3,Apicomplexa=28,Euglenozoa=11,Metazoa=6)
    assert (d[['C','B','N','M']]>0).all().all()
    assert int(d.N.sum())==2299 and int(d.C.sum())==1973
    rng=np.random.default_rng(20260922);draws={};raw={}
    for clade in CLADES:
        z=d[d.clade==clade].sort_values('species').reset_index(drop=True)
        a=[];b=[]
        for rep in range(1000):
            take=z.iloc[rng.choice(len(z),3,replace=False)]
            c,bg,n,m=take[['C','B','N','M']].sum()
            a.append((c/bg,(n/c)/(m/bg),n/m))
            b.append((c/3,n/c,n/3))
        draws[clade]=np.array(a);raw[clade]=np.array(b)
    rows=[]
    for scope,store in [('background_normalized',draws),('unadjusted',raw)]:
        for clade in CLADES[:-1]:
            ratios=store[clade]/store['Metazoa']
            np.testing.assert_allclose(ratios[:,0]*ratios[:,1],ratios[:,2])
            for k,component in enumerate(COMPONENTS):
                x=ratios[:,k]
                rows.append(dict(scope=scope,clade=clade,component=component,
                    geometric_mean=float(np.exp(np.log(x).mean())),q025=float(np.quantile(x,.025)),
                    q975=float(np.quantile(x,.975)),draws=1000,species_per_draw=3))
    result=pd.DataFrame(rows)
    expected=pd.read_csv(HERE/'expected/normalized_contrasts.tsv',sep='\t')
    chk=result[result.scope=='background_normalized'].merge(expected,on=['clade','component'],suffixes=('_computed','_expected'),validate='one_to_one')
    assert len(chk)==9
    maxerr=0.
    for col in ['geometric_mean','q025','q975']:
        np.testing.assert_allclose(chk[col+'_computed'],chk[col+'_expected'],rtol=1e-11,atol=1e-12)
        maxerr=max(maxerr,float(abs(chk[col+'_computed']-chk[col+'_expected']).max()))
    result.to_csv(out/'quantity_comparison.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(8,3.7),layout='constrained')
    for ax,scope,title in zip(axes,['unadjusted','background_normalized'],['Unadjusted counts','Normalized to genomic background']):
        for j,(comp,label,color) in enumerate(zip(COMPONENTS[:2],['HGT gene families','Genes per HGT family'],['#4477AA','#CC6677'])):
            z=result[(result.scope==scope)&(result.component==comp)].set_index('clade').loc[CLADES[:-1]]
            y=z.geometric_mean.to_numpy();ax.errorbar(np.arange(3)+(j-.5)*.16,y,yerr=[y-z.q025.to_numpy(),z.q975.to_numpy()-y],fmt='o',capsize=3,color=color,label=label)
        ax.set_yscale('log',base=2);ax.axhline(1,color='.6',ls=':',lw=.8)
        ax.set_xticks(range(3),['Amoebozoa','Apicomplexa','Euglenozoa'],rotation=15)
        ax.set_title(title);ax.set_ylabel('Fold difference from sampled nematodes')
        ax.spines[['top','right']].set_visible(False)
    axes[0].legend(frameon=False,fontsize=8)
    for ext in ['pdf','svg']:fig.savefig(out/('quantity_example.'+ext))
    plt.close(fig)
    validation=dict(status='passed',species=48,candidate_loci=2299,species_candidate_family_units=1973,
        saved_contrasts_verified=9,max_absolute_error=maxerr,seed=20260922,replicates=1000,
        interval_definition='central 95% species-subset range; not a confidence interval or formal P value',
        python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__)
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation,indent=2))
if __name__=='__main__':main()
