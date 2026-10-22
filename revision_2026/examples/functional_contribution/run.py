#!/usr/bin/env python3
"""Reproduce functional contributions from frozen matched-family profiles.
Matching and homology assignment are upstream inputs, not rerun here.
"""
from pathlib import Path
import argparse,json,platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;CATS=list('EFGH')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=HERE/'output')
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    a=pd.read_csv(HERE/'input/family_effects.tsv.gz',sep='\t')
    assert len(a)==1940 and a.species.nunique()==48 and not a.duplicated(['species','OG']).any()
    nf=a.groupby('species').size();ns=a.groupby('clade').species.nunique()
    weight=1/a.species.map(nf)/a.clade.map(ns)/4
    np.testing.assert_allclose(weight,a.weight);a['weight']=weight
    assert np.isclose(weight.sum(),1)
    summaries=[];curves=[];deletions=[];props=[]
    for cat in CATS:
        d=a['candidate_'+cat]-a['background_'+cat]
        np.testing.assert_allclose(d,a[cat],atol=1e-14)
        con=pd.DataFrame({'OG':a.OG,'signed':weight*d,'positive':weight*d.clip(lower=0),'negative':weight*d.clip(upper=0)})
        by=con.groupby('OG',as_index=False)[['signed','positive','negative']].sum()
        positive=by[by.positive>0].sort_values(['positive','OG'],ascending=[False,True]).copy()
        positive['rank']=np.arange(1,len(positive)+1);positive['cumulative_positive_share']=positive.positive.cumsum()/positive.positive.sum();positive['category']=cat
        curves.append(positive)
        gross=positive.positive.sum();net=by.signed.sum()
        summaries.append(dict(category=cat,positive_OGs=len(positive),n50=int(positive.loc[positive.cumulative_positive_share>=.5,'rank'].iloc[0]),n80=int(positive.loc[positive.cumulative_positive_share>=.8,'rank'].iloc[0]),top5_share=float(positive.positive.head(5).sum()/gross),gross_positive_pp=gross*100,gross_negative_pp=by.negative.sum()*100,net_effect_pp=net*100))
        props.append(dict(category=cat,candidate_percent=float((weight*a['candidate_'+cat]).sum()*100),background_percent=float((weight*a['background_'+cat]).sum()*100),difference_pp=net*100))
        largest=a.groupby('OG').loci.sum().reset_index().sort_values(['loci','OG'],ascending=[False,True])
        rules={'remove_top5_locus_OGs':set(largest.head(5).OG),'remove_top5_category_positive_OGs':set(positive.head(5).OG)}
        baseline=a.groupby(['species','clade'],as_index=False)[cat].mean()
        for rule,drop in rules.items():
            rem=a[~a.OG.isin(drop)].groupby(['species','clade']).agg(n=('OG','size'),effect=(cat,'mean')).reset_index()
            rem=rem[rem.n>=10];counts=rem.groupby('clade').species.nunique();rem=rem[rem.clade.isin(counts[counts>=3].index)]
            base=baseline[baseline.species.isin(rem.species)]
            deletions.append(dict(category=cat,deletion=rule,retained_species=rem.species.nunique(),retained_clades=rem.clade.nunique(),baseline_same_species_pp=base.groupby('clade')[cat].mean().mean()*100,remaining_effect_pp=rem.groupby('clade').effect.mean().mean()*100,baseline_all_species_pp=net*100))
    su=pd.DataFrame(summaries);de=pd.DataFrame(deletions);cu=pd.concat(curves,ignore_index=True);pr=pd.DataFrame(props)
    errors=[]
    for got,fn,keys,cols in [(su,'contribution_summary.tsv',['category'],['positive_OGs','n50','n80','top5_share','gross_positive_pp','gross_negative_pp','net_effect_pp']),(de,'deletion_summary.tsv',['category','deletion'],['retained_species','retained_clades','baseline_same_species_pp','remaining_effect_pp','baseline_all_species_pp'])]:
        exp=pd.read_csv(HERE/'expected'/fn,sep='\t');j=got.merge(exp,on=keys,suffixes=('_computed','_expected'),validate='one_to_one');assert len(j)==len(exp)==len(got)
        for c in cols:
            np.testing.assert_allclose(j[c+'_computed'],j[c+'_expected'],atol=1e-10,rtol=1e-10)
            errors.append(float(abs(j[c+'_computed']-j[c+'_expected']).max()))
    for frame,fn in [(su,'contribution_summary.tsv'),(de,'deletion_summary.tsv'),(cu,'contribution_curves.tsv'),(pr,'functional_proportions.tsv')]:frame.to_csv(out/fn,sep='\t',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(6,3.8),layout='constrained')
    for cat,label,col in zip(CATS,['Amino-acid','Nucleotide','Carbohydrate','Coenzyme'],['#4477AA','#EE6677','#228833','#AA3377']):
        z=cu[cu.category==cat];ax.plot(np.r_[0,z['rank']],np.r_[0,z.cumulative_positive_share]*100,label=label,color=col)
    ax.axhline(50,c='.6',ls=':',lw=.8);ax.set_xlabel('HGT gene families ranked by positive contribution')
    ax.set_ylabel('Cumulative positive contribution (%)');ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
    for ext in ['pdf','svg']:fig.savefig(out/('functional_contribution_example.'+ext))
    plt.close(fig)
    validation=dict(status='passed',species=48,matched_species_family_observations=1940,categories_verified=4,deletion_comparisons_verified=8,max_absolute_error=max(errors),python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__,interpretation='Annotation contributions; not enzyme activity, pathway flux, or gene essentiality.')
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n');print(json.dumps(validation,indent=2))
if __name__=='__main__':main()
