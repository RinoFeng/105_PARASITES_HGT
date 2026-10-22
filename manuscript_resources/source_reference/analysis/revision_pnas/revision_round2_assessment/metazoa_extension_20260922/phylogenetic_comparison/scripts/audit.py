from pathlib import Path
import copy, hashlib, json
import numpy as np
import pandas as pd
from Bio import Phylo

OUT = Path(__file__).resolve().parents[1]
B = OUT.parent
R = B.parent
P = R/'innovation_execution/family_comparison/validation_followup/species_phylogeny'
CLADES = ['Apicomplexa', 'Euglenozoa', 'Amoebozoa']

def inputs():
    s=pd.read_csv(B/'statistical_feasibility/species_endpoints48.tsv',sep='\t')
    m=pd.read_csv(P/'species_manifest.tsv',sep='\t')
    s=s.merge(m[['species','species_id','single_complete']],on='species',validate='one_to_one')
    s=s.merge(pd.read_csv(P/'species_marker_coverage.tsv',sep='\t')[['species','retained_markers']],on='species',validate='one_to_one')
    return s

def limited_tree(frame, model='LG', fraction=.5):
    t=Phylo.read(P/f'species_{model}.treefile','newick')
    keep=set(frame.species_id)
    assert keep <= {x.name for x in t.get_terminals()}
    for x in list(t.get_terminals()):
        if x.name not in keep: t.prune(x)
    nem=list(frame.loc[frame.clade=='Metazoa','species_id'])
    # Place a computational root on the edge separating the two sampled clades.
    # This is NOT an inference about the root of eukaryotes.
    distances=np.array([[t.distance(a,b) for b in frame.species_id] for a in frame.species_id])
    nemtips=[x for x in t.get_terminals() if x.name in nem]
    if not t.is_monophyletic(nemtips):
        other=list(frame.loc[frame.clade!='Metazoa','species_id'])
        t.root_with_outgroup(other)
    assert t.is_monophyletic(nemtips)
    t.root_with_outgroup(nem, outgroup_branch_length=0)
    assert len(t.root.clades)==2
    total=sum(x.branch_length or 0 for x in t.root.clades)
    for c in t.root.clades:
        c.branch_length=total*(fraction if {x.name for x in c.get_terminals()}==set(nem) else 1-fraction)
    t.root.branch_length=0
    assert np.allclose(distances,[[t.distance(a,b) for b in frame.species_id] for a in frame.species_id])
    return t

def covariance(t, ids, scale='substitution', lam=1):
    t=copy.deepcopy(t)
    if scale=='grafen':
        heights={c:(len(c.get_terminals())-1)/(len(t.get_terminals())-1) for c in t.find_clades()}
        for c in t.find_clades():
            for child in c.clades: child.branch_length=heights[c]-heights[child]
        t.root.branch_length=0
    d=t.depths()
    leaves={x.name:x for x in t.get_terminals()}
    v=np.array([[d[t.common_ancestor(leaves[a],leaves[b])] for b in ids] for a in ids])
    v=lam*v+(1-lam)*np.diag(np.diag(v))
    v/=np.mean(np.diag(v))
    assert np.min(np.linalg.eigvalsh(v))>0
    return v

def main():
    OUT.mkdir(exist_ok=True)
    s=inputs()
    f=pd.read_csv(R/'family_contribution_extension_20260922/tables/family_profiles47.tsv.gz',sep='\t')
    a=pd.read_csv(B/'tables/recovered_family_profiles.tsv.gz',sep='\t')
    a=a[(a.species=='Angiostrongylus costaricensis') & a.same_original_reference]
    f=pd.concat([f,a[f.columns]],ignore_index=True)
    assert not f.duplicated(['species','OG','role']).any()
    f.to_csv(OUT/'family_inputs48.tsv.gz',sep='\t',index=False)
    metrics=['eggnog_fraction','informative_cog_fraction','ko_fraction','pfam_fraction','length']
    rows=[]
    for (sp,role),g in f.groupby(['species','role']):
        row={'species':sp,'role':role,'families':len(g),'loci':g.units.sum()}
        row.update({x:g[x].mean() for x in metrics}); rows.append(row)
    pool=pd.DataFrame(rows);pool.to_csv(OUT/'annotation_by_pool.tsv',sep='\t',index=False)
    for role in ['candidate','background']:
        z=pool[pool.role==role].drop(columns='role').rename(columns={x:role+'_'+x for x in ['families','loci']+metrics})
        s=s.merge(z,on='species',validate='one_to_one')
    assert (s.C==s.candidate_families).all() and (s.B==s.background_families).all()
    assert (s.N==s.candidate_loci).all() and (s.M==s.background_loci).all()
    ledger=pd.read_csv(R/'innovation_final8989_20260915/macro/record_locus_family_quarantine_ledger.tsv',sep='\t')
    ledger=ledger[ledger.final_catalogue].copy()
    ledger['strict']=ledger.direct_evidence_status.eq('own_sequence_exact_eligible_representative')
    recovery=pd.read_csv(B/'tables/candidate_locus_audit.tsv',sep='\t')
    rec=recovery[recovery.same_original_reference & recovery.in_final_catalogue].set_index('record_id')
    ix=ledger.record_id.isin(rec.index)
    ledger.loc[ix,'strict']=ledger.loc[ix,'record_id'].map(rec.status.eq('own_sequence_exact_eligible_representative')).values
    q=ledger.groupby('species').agg(catalogue_records=('record_id','size'),strict_records=('strict','sum')).reset_index()
    q['candidate_record_link_rate']=q.strict_records/q.catalogue_records
    s=s.merge(q,on='species',validate='one_to_one')
    # Link rate is outcome/selection related: diagnostic, never a routine adjustment covariate.
    s.to_csv(OUT/'species_observation_audit.tsv',sep='\t',index=False)
    trees=[];splitrows=[];overlap=[]
    for clade in CLADES:
        z=s[s.clade.isin([clade,'Metazoa'])].copy()
        for model in ['LG','PMSF']:
            t=limited_tree(z,model)
            Phylo.write(t,OUT/f'{clade}_{model}_limited.nwk','newick')
            v=covariance(t,list(z.species_id))
            trees.append(dict(comparison=clade,model=model,species=len(z),focal_genera=z[z.clade==clade].genus.nunique(),nematode_genera=6,all_tips_present=True,clade_split_present=True,positive_branches=all((c.branch_length or 0)>0 for c in t.find_clades() if c!=t.root),condition_number=np.linalg.cond(v),min_markers=z.retained_markers.min(),max_markers=z.retained_markers.max(),branch_scale='substitutions_per_site_not_time'))
            for genus,g in z.groupby('genus'):
                ids=list(g.species_id)
                splitrows.append(dict(comparison=clade,model=model,genus=genus,species=len(ids),monophyletic=bool(t.is_monophyletic([x for x in t.get_terminals() if x.name in ids]))))
        for metric in ['retained_markers','single_complete','background_eggnog_fraction','background_ko_fraction','background_informative_cog_fraction','background_length','candidate_record_link_rate']:
            x=z.loc[z.clade==clade,metric]; y=z.loc[z.clade=='Metazoa',metric]
            lo=max(x.min(),y.min()); hi=min(x.max(),y.max())
            overlap.append(dict(comparison=clade,covariate=metric,focal_min=x.min(),focal_max=x.max(),nematode_min=y.min(),nematode_max=y.max(),shared_min=lo,shared_max=hi,overlap_exists=lo<=hi,focal_in_range=x.between(lo,hi).sum(),nematode_in_range=y.between(lo,hi).sum()))
    pd.DataFrame(trees).to_csv(OUT/'tree_audit.tsv',sep='\t',index=False)
    pd.DataFrame(splitrows).to_csv(OUT/'genus_tree_audit.tsv',sep='\t',index=False)
    pd.DataFrame(overlap).to_csv(OUT/'observation_overlap.tsv',sep='\t',index=False)
    print(pd.DataFrame(trees).to_string(index=False))
    print(pd.DataFrame(overlap).to_string(index=False))
    src=[B/'statistical_feasibility/EXPLORATORY_TEST_PLAN.json',B/'statistical_feasibility/species_endpoints48.tsv',P/'species_LG.treefile',P/'species_PMSF.treefile',R/'family_contribution_extension_20260922/tables/family_profiles47.tsv.gz',B/'tables/recovered_family_profiles.tsv.gz']
    pd.DataFrame([{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in src]).to_csv(OUT/'input_manifest.tsv',sep='\t',index=False)

if __name__=='__main__': main()
