"""Family-level decomposition of existing broad functional annotations.

No sequence analysis, origin relabeling, mechanistic prioritization or causal test.
"""
import csv,gzip,json,re,hashlib,itertools
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
import pandas as pd

B=Path(__file__).resolve().parent; A=Path('PROJECT_ROOT/analysis/revision_pnas/revision_round2_assessment/innovation_execution/source_function_convergence_20260914'); O=Path('PROJECT_ROOT/analysis/revision_pnas/revision_round2_assessment/innovation_execution'); G=Path('PROJECT_ROOT/analysis/revision_pnas/revision_round2_assessment/innovation_execution/genome_integration'); R=G/'refinement_20260914'; OLD=B.parent/'function_robustness'
CATS=list('EFGH'); PRIMARY=['Amoebozoa','Apicomplexa','Euglenozoa','Metazoa'];SEED=20260914
CHECKS=[];INPUTS={}
def rows(p):
    INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as h:
        yield from csv.DictReader(h,delimiter='\t')
def write(name,rr):
    if isinstance(rr,pd.DataFrame):df=rr
    else:df=pd.DataFrame(rr)
    df.to_csv(B/name,sep='\t',index=False,compression='infer')
def check(name,b):
    CHECKS.append((name,bool(b)))
    if not b:raise AssertionError(name)
def truth(x):return str(x).lower() in ['true','1','yes']
def toks(x):return [z for z in re.split(r'[;|,]',x or '') if z and z not in ['NA','-']]
def labels(r):return sorted(set(r.get('COG_category',''))&set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')-set('RS'))
def meanstats(vals):
    return dict(mean=float(np.mean(vals)),median=float(np.median(vals)),minimum=float(np.min(vals)),maximum=float(np.max(vals)))

def main():
    freeze=json.loads((B/'stage_freeze.json').read_text())
    check('stage_groups_frozen',freeze['membership_sha256']==hashlib.sha256((B/'stage_case_membership_frozen.tsv').read_bytes()).hexdigest())
    eco={r['species']:r for r in rows(A/'ecology_remaining_review_20260914/eligible_species_ecology_review2_52.tsv')}
    ss=set(eco); check('52_species',len(ss)==52)
    pairs=defaultdict(list)
    for r in rows(OLD/'annotation_matched_pairs.tsv.gz'):
        if r['species'] in ss and r['mode']=='length_fine_annotation':pairs[(r['species'],r['candidate_OG'])].append(r['background_OG'])
    needed={(s,o,'candidate') for s,o in pairs}|{(s,b,'background') for (s,o),bs in pairs.items() for b in bs}
    metrics={}
    for r in rows(OLD/'family_annotation_metrics.tsv.gz'):
        key=(r['species'],r['OG'],r['role'])
        if key in needed:metrics[key]=r
    check('all_needed_family_metrics',len(metrics)==len(needed))
    strict={r['unit_id']:r for r in rows(B.parent/'structure/strict_candidate_features.tsv')}
    og={r['unit_id']:r['inferred_family'] for r in rows(R/'unit_orthogroup_membership.tsv.gz')}
    blocked={r['inferred_family'] for r in rows(R/'resolved_family_function_audit.tsv.gz') if truth(r['excluded_ordinary_background'])}
    possible=set()
    for r in rows(R/'candidate_locus_evidence_audit.tsv'):
        possible.update(toks(r['own_record_gene_ids'])+toks(r['coordinate_unit_ids'])+toks(r['legacy_unit_ids']))
    count=Counter();sums=defaultdict(lambda:np.zeros(4));kos=defaultdict(lambda:[set() for _ in CATS]);quality=defaultdict(Counter);treeclasses=defaultdict(Counter)
    for inp in rows(G/'species_inputs.tsv'):
        sp=inp['species']
        if sp not in ss:continue
        for r in rows(G/'measurements'/(inp['species_id']+'.tsv.gz')):
            uid=r['unit_id'];family=og[uid]
            role='candidate' if uid in strict else ('background' if truth(r['eligible_structure']) and r['analysis_role']=='non_catalogue_background' and uid not in possible and family not in blocked else None)
            key=(sp,family,role)
            if key not in needed:continue
            count[key]+=1;ll=labels(r)
            if ll:
                for i,cat in enumerate(CATS):
                    if cat in ll:
                        sums[key][i]+=1/len(ll)
                        kos[key][i].update(re.findall(r'K\d{5}',r.get('KEGG_ko','')))
            if role=='candidate':
                z=strict[uid];quality[sp]['matched_candidate_loci']+=1
                quality[sp]['any_Tier1_label']+=('Tier1' in z['all_contig_tiers'])
                quality[sp]['legacy_explicit_no_long_branch']+=(z['all_legacy_long_branch'].strip().upper()=='FALSE')
                treeclasses[sp][z['original_tree_classes']]+=1
        print('ANNOTATION_SUMMARY_READ',sp,flush=True)
    profiles={k:sums[k]/count[k] for k in needed}
    for k,m in metrics.items():check('member_count_'+str(k),count[k]==int(m['units']))
    family_rows=[];objects=[]
    species_fam=Counter(s for s,o in pairs)
    for (s,o),bs in sorted(pairs.items()):
        check('valid_control_count_'+s+o,3<=len(bs)<=5 and len(set(bs))==len(bs))
        p=profiles[(s,o,'candidate')];q=np.stack([profiles[(s,b,'background')] for b in bs]);delta=p-q.mean(axis=0)
        objects.append(dict(species=s,OG=o,clade=eco[s]['prespecified_clade'],state=eco[s]['review2_cell_state'],P=p,Q=q,D=delta,controls=bs,units=count[(s,o,'candidate')]))
        for i,c in enumerate(CATS):family_rows.append(dict(species=s,clade=eco[s]['prespecified_clade'],state=eco[s]['review2_cell_state'],OG=o,category=c,
            candidate_share=p[i],background_share=q[:,i].mean(),difference=delta[i],candidate_loci=count[(s,o,'candidate')],matched_background_families=len(bs)))
    df=pd.DataFrame(family_rows)
    write('matched_family_function_effects.tsv.gz',df)
    species=df.groupby(['species','clade','state','category'],as_index=False).agg(candidate_share=('candidate_share','mean'),background_share=('background_share','mean'),effect=('difference','mean'),families=('OG','nunique'))
    old={(r['species'],r['category']):r for r in rows(OLD/'direction_species_effects.tsv') if r['scenario']=='all_members_length_fine_annotation' and r['species'] in ss}
    check('208_species_function_rows',len(species)==208)
    for r in species.to_dict('records'):
        prior=old[(r['species'],r['category'])]
        for key in ['candidate_share','background_share','effect']:check('reproduce_'+r['species']+r['category']+key,abs(r[key]-float(prior[key]))<1e-12)
        check('family_number_'+r['species']+r['category'],r['families']==int(prior['families']))
    write('species_function_atlas_data.tsv',species)
    reliability=[]
    for s in sorted(ss):
        mm=[metrics[(sp,o,'candidate')] for sp,o in pairs if sp==s];qq=quality[s];n=qq['matched_candidate_loci']
        reliability.append(dict(species=s,clade=eco[s]['prespecified_clade'],state=eco[s]['review2_cell_state'],matched_families=len(mm),matched_candidate_loci=n,
            informative_COG_fraction_family_mean=np.mean([float(m['informative_cog_fraction']) for m in mm]),
            stronger_annotation_fraction_family_mean=np.mean([float(m['stronger_match_fraction']) for m in mm]),
            any_Tier1_label_locus_fraction=qq['any_Tier1_label']/n,legacy_explicit_no_long_branch_locus_fraction=qq['legacy_explicit_no_long_branch']/n,
            scope='Matched strict-CDS subset; these are separate evidence indicators, not a calibrated reliability grade.'))
    write('sample_reliability_indicators.tsv',reliability)
    write('preserved_original_tree_class_counts.tsv',[dict(species=s,original_tree_classes=k,matched_loci=v) for s,cc in treeclasses.items() for k,v in cc.items()])
    ogspecies=defaultdict(set);ogclades=defaultdict(set);ogloci=Counter()
    for obj in objects:
        ogspecies[obj['OG']].add(obj['species']);ogclades[obj['OG']].add(obj['clade']);ogloci[obj['OG']]+=obj['units']
    sharing={o:('cross_clade_shared' if len(ogclades[o])>1 else 'within_clade_multispecies' if len(ogspecies[o])>1 else 'single_species_observed') for o in ogspecies}
    write('OG_observation_scope.tsv',[dict(OG=o,species=len(ogspecies[o]),clades=len(ogclades[o]),matched_candidate_loci=ogloci[o],observation_class=sharing[o],
         clade_labels=';'.join(sorted(ogclades[o])),independent_transfers_confirmed=False) for o in sorted(ogspecies)])
    ns=Counter(eco[s]['prespecified_clade'] for s in ss)
    contrib=[]
    for r in family_rows:
        cl=r['clade'];w=1/species_fam[r['species']]/ns[cl]
        contrib.append(dict(r,observation_class=sharing[r['OG']],clade_weight=w,clade_contribution=w*r['difference'],
              primary_weight=w/4 if cl in PRIMARY else 0,primary_contribution=w*r['difference']/4 if cl in PRIMARY else 0,
              primary_positive_part=w*max(r['difference'],0)/4 if cl in PRIMARY else 0,
              primary_negative_part=w*min(r['difference'],0)/4 if cl in PRIMARY else 0))
    cd=pd.DataFrame(contrib);write('family_contribution_decomposition.tsv.gz',cd)
    byog=cd.groupby(['category','OG','observation_class'],as_index=False)[['primary_contribution','primary_positive_part','primary_negative_part']].sum()
    byog['matched_loci']=byog.OG.map(ogloci);byog['observed_species']=byog.OG.map(lambda o:len(ogspecies[o]));byog['observed_clades']=byog.OG.map(lambda o:len(ogclades[o]))
    write('OG_function_contributions.tsv',byog)
    decomp=cd.groupby(['clade','category','observation_class'],as_index=False).agg(contribution=('clade_contribution','sum'),OGs=('OG','nunique'))
    write('clade_sharing_decomposition.tsv',decomp)
    primary=byog.groupby(['category','observation_class'],as_index=False)[['primary_contribution','primary_positive_part','primary_negative_part']].sum();write('primary_sharing_decomposition.tsv',primary)
    top_loci=set(sorted(ogloci,key=lambda o:(-ogloci[o],o))[:5]);top_cat={}
    dominance=[]
    for c in CATS:
        z=byog[byog.category==c].sort_values(['primary_positive_part','OG'],ascending=[False,True]);top_cat[c]=set(z.head(5).OG)
        pos=z.primary_positive_part.sum()
        dominance.append(dict(category=c,top5_positive_contribution_share=z.head(5).primary_positive_part.sum()/pos if pos else np.nan,
            top1_positive_contribution_share=z.head(1).primary_positive_part.sum()/pos if pos else np.nan,
            gross_positive=pos,gross_negative=z.primary_negative_part.sum(),net_effect=z.primary_contribution.sum(),
            top5_OGs=';'.join(z.head(5).OG),top5_locus_OGs=';'.join(sorted(top_loci))))
        base=species[(species.category==c)&species.clade.isin(PRIMARY)].groupby('clade').effect.mean().mean()
        check('exact_primary_decomposition_'+c,abs(base-z.primary_contribution.sum())<1e-12)
    write('dominance_summary.tsv',dominance)
    deletion=[];deletion_species=[]
    for c in CATS:
        z=df[df.category==c];base=species[species.category==c].set_index('species')
        for name,remove in [('remove_top5_locus_OGs',top_loci),('remove_top5_category_positive_OGs',top_cat[c]),('exclude_cross_clade_shared_OGs',{o for o in sharing if sharing[o]=='cross_clade_shared'})]:
            zz=z[~z.OG.isin(remove)];gg=zz.groupby('species').agg(n=('OG','nunique'),effect=('difference','mean'))
            gg=gg[gg.n>=10];gg['clade']=[eco[s]['prespecified_clade'] for s in gg.index]
            common=[cl for cl in PRIMARY if sum(gg.clade==cl)>=3]
            for s,r in gg.iterrows():deletion_species.append(dict(scenario=name,category=c,species=s,clade=r.clade,families=r.n,baseline_effect=base.loc[s,'effect'],remaining_effect=r.effect))
            for cl in PRIMARY:
                gs=gg[gg.clade==cl];bn=base[base.clade==cl];valid=len(gs)>=3
                deletion.append(dict(scenario=name,category=c,clade=cl,remaining_species=len(gs),baseline_species=len(bn),comparable=valid,
                    baseline_same_species=base.loc[gs.index,'effect'].mean() if len(gs) else np.nan,
                    remaining_effect=gs.effect.mean() if len(gs) else np.nan,removed_OGs=len(remove)))
            keep=gg[gg.clade.isin(common)]
            deletion.append(dict(scenario=name,category=c,clade='PRIMARY_EQUAL_RETAINED_CLADES',remaining_species=len(keep),baseline_species=sum(ns[x] for x in PRIMARY),
                comparable=len(common)>=3,baseline_same_species=base.loc[keep.index].groupby('clade').effect.mean().mean() if len(keep) else np.nan,
                remaining_effect=keep.groupby('clade').effect.mean().mean() if len(keep) else np.nan,removed_OGs=len(remove)))
    write('deletion_same_fine_background.tsv',deletion);write('deletion_species_effects.tsv',deletion_species)
    # Fixed matched controls: global OG ranks keep repeated backgrounds coupled.
    rng=np.random.default_rng(SEED);bgnames=sorted({b for o in objects for b in o['controls']});bi={b:i for i,b in enumerate(bgnames)}
    ranks=rng.random((500,len(bgnames)));draw_species=defaultdict(list)
    for obj in objects:
        ids=[bi[b] for b in obj['controls']];chosen=ranks[:,ids].argmin(axis=1)
        draw_species[obj['species']].append(obj['P'][None,:]-obj['Q'][chosen])
    draw_species={s:np.mean(v,axis=0) for s,v in draw_species.items()}
    background_sensitivity=[]
    for cl in PRIMARY:
        draw=np.mean([v for s,v in draw_species.items() if eco[s]['prespecified_clade']==cl],axis=0)
        for i,c in enumerate(CATS):background_sensitivity.append(dict(clade=cl,category=c,draw_q05=np.quantile(draw[:,i],.05),draw_q95=np.quantile(draw[:,i],.95),
             positive_draw_fraction=np.mean(draw[:,i]>0),replicates=500,meaning='Sensitivity to choosing one of fixed matched controls, not a significance test or full confidence interval.'))
    write('matched_background_choice_sensitivity.tsv',background_sensitivity)
    # Category-linked OG and KO overlap. All clade comparisons retained.
    aset=defaultdict(set);ak=defaultdict(set);qset=defaultdict(set);qk=defaultdict(set)
    for obj in objects:
        cl=obj['clade'];s=obj['species'];o=obj['OG']
        for i,c in enumerate(CATS):
            if obj['P'][i]>0:aset[(cl,c)].add(o);ak[(cl,c)].update(kos[(s,o,'candidate')][i])
            for b in obj['controls']:
                if profiles[(s,b,'background')][i]>0:qset[(cl,c)].add(b);qk[(cl,c)].update(kos[(s,b,'background')][i])
    overlap=[]
    def jac(x,y):return len(x&y)/len(x|y) if x|y else np.nan
    for x,y in itertools.combinations(PRIMARY,2):
        for c in CATS:
            u,v=aset[(x,c)],aset[(y,c)];ku,kv=ak[(x,c)],ak[(y,c)]
            bg1,bg2=qset[(x,c)],qset[(y,c)];bk1,bk2=qk[(x,c)],qk[(y,c)]
            shared=u&v;xu=u-shared;yv=v-shared
            # Shared KO attached to nonoverlapping OGs, only category-contributing genes.
            ku_dist=set();kv_dist=set()
            for ob in objects:
                if ob['clade']==x and ob['OG'] in xu:ku_dist.update(kos[(ob['species'],ob['OG'],'candidate')][CATS.index(c)])
                if ob['clade']==y and ob['OG'] in yv:kv_dist.update(kos[(ob['species'],ob['OG'],'candidate')][CATS.index(c)])
            overlap.append(dict(clade1=x,clade2=y,category=c,candidate_OGs_1=len(u),candidate_OGs_2=len(v),shared_OGs=len(shared),OG_jaccard=jac(u,v),
                 candidate_KOs_1=len(ku),candidate_KOs_2=len(kv),shared_KOs=len(ku&kv),KO_jaccard=jac(ku,kv),
                 shared_KOs_on_nonoverlapping_OGs=len(ku_dist&kv_dist),
                 matched_background_OG_jaccard=jac(bg1,bg2),matched_background_KO_jaccard=jac(bk1,bk2),
                 interpretation='Set overlaps; background union has different set sizes, no calibrated enrichment or independent convergence claim.'))
    write('cross_clade_OG_KO_overlap.tsv',overlap)
    # Same number of candidate-family slots; sample matched backgrounds to assess overlap scale.
    baseline_pairs={(r['clade1'],r['clade2'],r['category']):r for r in overlap}
    sim=defaultdict(list); ko_sim=defaultdict(list)
    for rep in range(200):
        oo=defaultdict(set);kk=defaultdict(set)
        for ob in objects:
            b=ob['controls'][int(ranks[rep,[bi[b] for b in ob['controls']]].argmin())];key=(ob['species'],b,'background')
            for i,c in enumerate(CATS):
                if profiles[key][i]>0:oo[(ob['clade'],c)].add(b);kk[(ob['clade'],c)].update(kos[key][i])
        for x,y in itertools.combinations(PRIMARY,2):
            for c in CATS:
                sim[(x,y,c)].append(jac(oo[(x,c)],oo[(y,c)]));ko_sim[(x,y,c)].append(jac(kk[(x,c)],kk[(y,c)]))
    simrows=[]
    for key,values in sim.items():
        x,y,c=key
        for label,vv in [('OG',values),('KO',ko_sim[key])]:
            vv=np.asarray(vv);vv=vv[np.isfinite(vv)]
            simrows.append(dict(clade1=x,clade2=y,category=c,metric=label,candidate_jaccard=baseline_pairs[key][label+'_jaccard'],
                background_draw_median=np.median(vv) if len(vv) else np.nan,background_q05=np.quantile(vv,.05) if len(vv) else np.nan,
                background_q95=np.quantile(vv,.95) if len(vv) else np.nan,valid_draws=len(vv),draws=200,
                limitation='Equal number of matched family slots, not equal KO-bearing family count or family member count; descriptive background benchmark only.'))
    write('overlap_background_draw_benchmark.tsv',simrows)
    # Audit additional member-count comparability without changing primary matches.
    sizeaudit=[]
    for ob in objects:
        n=ob['units'];good=[b for b in ob['controls'] if .5<=count[(ob['species'],b,'background')]/n<=2]
        sizeaudit.append(dict(species=ob['species'],OG=ob['OG'],candidate_loci=n,original_controls=len(ob['controls']),member_count_ratio_half_to_two_controls=len(good),at_least_three=len(good)>=3))
    write('family_member_count_matching_audit.tsv',sizeaudit)
    sizeobjects=[]
    for ob in objects:
        bs=[b for b in ob['controls'] if .5<=count[(ob['species'],b,'background')]/ob['units']<=2]
        if len(bs)>=3:
            sizeobjects.append(dict(ob,Q=np.stack([profiles[(ob['species'],b,'background')] for b in bs])))
    size_species=Counter(ob['species'] for ob in sizeobjects)
    sizegroups=defaultdict(list)
    for ob in sizeobjects:
        if size_species[ob['species']]>=10:sizegroups[ob['species']].append(ob['P']-ob['Q'].mean(axis=0))
    size_rows=[]
    for s,vals in sizegroups.items():
        for i,c in enumerate(CATS):
            prior=old[(s,c)]
            size_rows.append(dict(species=s,clade=eco[s]['prespecified_clade'],category=c,families=size_species[s],
                baseline_effect=float(prior['effect']),member_count_matched_effect=np.mean(vals,axis=0)[i]))
    write('member_count_matched_species_sensitivity.tsv',size_rows)
    size_df=pd.DataFrame(size_rows)
    size_clades=size_df.groupby(['clade','category'],as_index=False).agg(species=('species','nunique'),baseline_same_species=('baseline_effect','mean'),member_count_matched_effect=('member_count_matched_effect','mean'))
    write('member_count_matched_clade_sensitivity.tsv',size_clades)
    # Functional summaries only after the stage membership file was frozen.
    membership=list(rows(B/'stage_case_membership_frozen.tsv'));members=pd.DataFrame(membership)
    stage=members.merge(species,on=['species','clade'],validate='one_to_many');write('stage_case_species_function.tsv',stage)
    caserows=[]
    for (case,group,c),t in stage.groupby(['case','group','category']):
        vals=t.effect.to_numpy();loo=[np.delete(vals,i).mean() for i in range(len(vals))] if len(vals)>1 else []
        caserows.append(dict(case=case,group=group,category=c,species=len(t),mean_effect=np.mean(vals),median_effect=np.median(vals),
            positive_species=sum(vals>1e-12),leave_one_species_min=min(loo) if loo else np.nan,leave_one_species_max=max(loo) if loo else np.nan,
            inference='Case description only; groups are not independent ecological transitions.'))
    write('stage_case_function_summary.tsv',caserows)
    contrast=[]
    for case,t in stage.groupby('case'):
        for g1,g2 in itertools.combinations(sorted(t.group.unique()),2):
            for c in CATS:
                a=t[(t.group==g1)&(t.category==c)];b=t[(t.group==g2)&(t.category==c)]
                contrast.append(dict(case=case,group1=g1,group2=g2,category=c,n1=len(a),n2=len(b),difference_2_minus_1=b.effect.mean()-a.effect.mean(),
                    both_at_least3=min(len(a),len(b))>=3,genus_confounding_retained=True,independent_ecological_effect_claim=False))
    write('stage_case_pairwise_descriptions.tsv',contrast)
    for key,h in INPUTS.items():check('input_unchanged_'+key,hashlib.sha256(Path(key).read_bytes()).hexdigest()==h)
    check('stage_groups_still_frozen',freeze['membership_sha256']==hashlib.sha256((B/'stage_case_membership_frozen.tsv').read_bytes()).hexdigest())
    write('input_manifest.tsv',[dict(path=p,sha256=h) for p,h in INPUTS.items()])
    write('safety_skipped_tasks.tsv',pd.DataFrame(columns=['task','reason','status']))
    (B/'analysis_status.json').write_text(json.dumps(dict(status='complete',species=52,matched_candidate_family_units=len(objects),unique_candidate_OGs=len(ogspecies),
        matched_candidate_loci=sum(ob['units'] for ob in objects),primary_clades=PRIMARY,checks=len(CHECKS),all_checks_passed=all(x[1] for x in CHECKS),
        original_inputs_unchanged=True,causal_models_fitted=False,independent_transfer_claims=False,safety_skipped=[]),indent=2)+'\n')
    (B/'validation_details.json').write_text(json.dumps(CHECKS)+'\n')
    print((B/'analysis_status.json').read_text());print(pd.DataFrame(dominance).to_string(index=False))
if __name__=='__main__':main()
