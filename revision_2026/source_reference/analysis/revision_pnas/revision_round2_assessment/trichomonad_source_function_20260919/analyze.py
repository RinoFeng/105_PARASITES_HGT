"""Reproduce and standardize two existing species-level annotation summaries.

Read-only inputs, frozen candidate membership, no sequence/tree inference.
"""
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

B = Path(__file__).resolve().parent
R = B.parent
F = R / 'innovation_final8989_20260915'
G = R / 'innovation_execution/genome_integration'
TAX = Path('EXTERNAL_DATA_ROOT/Data/DB/taxonomy/tax_20230804')
SP = ['Trichomonas vaginalis G3', 'Tritrichomonas foetus']
SHORT = ['T. vaginalis', 'T. foetus']
CATS = list('EFGH')
STAGES = ['original', 'common', 'standardized']
INPUTS = set()
CHECKS = []
SEED = 20260919


def rows(p):
    p = Path(p)
    INPUTS.add(p)
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as h:
        yield from csv.DictReader(h, delimiter='\t')


def write(name, records, fields=None):
    df = records if isinstance(records, pd.DataFrame) else pd.DataFrame(records, columns=fields)
    df.to_csv(B / name, sep='\t', index=False, na_rep='NA')


def check(name, condition):
    CHECKS.append({'check': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def truth(x):
    return str(x).lower() in {'true', '1', 'yes'}


def splitids(x):
    return [s for s in (x or '').split(';') if s and s not in {'NA', '-'}]


def integer(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def read_profiles():
    original = {}
    for r in rows(F / 'ecology/matched_family_function_effects.tsv.gz'):
        if r['species'] not in SP:
            continue
        key = (r['species'], r['OG'])
        z = original.setdefault(key, {'P': np.zeros(4), 'B': np.zeros(4), 'units': int(r['candidate_loci'])})
        j = CATS.index(r['category'])
        z['P'][j] = float(r['candidate_share'])
        z['B'][j] = float(r['background_share'])
    check('original_family_counts', [sum(k[0] == s for k in original) for s in SP] == [321, 263])
    pairs = defaultdict(list)
    for r in rows(F / 'function_robustness/annotation_matched_pairs.tsv.gz'):
        if (r['species'], r['candidate_OG']) in original and r['mode'] in {'length_fine_annotation', 'length_only'}:
            pairs[(r['species'], r['candidate_OG'], r['mode'])].append(r['background_OG'])
    strict = {r['unit_id']: r for r in rows(F / 'structure/strict_candidate_features.tsv') if r['species'] in SP}
    blocked = {r['inferred_family'] for r in rows(G / 'refinement_20260914/resolved_family_function_audit.tsv.gz') if truth(r['excluded_ordinary_background'])}
    potential = set()
    for r in rows(G / 'refinement_20260914/candidate_locus_evidence_audit.tsv'):
        for field in ['own_record_gene_ids', 'coordinate_unit_ids', 'legacy_unit_ids']:
            potential.update(splitids(r[field]))
    measurements = []
    for code, species in [('S096', SP[0]), ('S098', SP[1])]:
        for r in rows(G / ('measurements/' + code + '.tsv.gz')):
            if truth(r['eligible_structure']):
                measurements.append((species, r['unit_id'], r['COG_category'], r['analysis_role']))
    wanted_units = {u for _, u, _, _ in measurements}
    unit_og = {r['unit_id']: r['inferred_family'] for r in rows(G / 'refinement_20260914/unit_orthogroup_membership.tsv.gz') if r['unit_id'] in wanted_units}
    needed = {(s, o, 'candidate') for s, o in original}
    needed.update((s, o, 'background') for (s, _, _), bs in pairs.items() for o in bs)
    sums = defaultdict(lambda: np.zeros(4))
    counts = Counter()
    family_units = defaultdict(list)
    for species, u, cog, role in measurements:
        og = unit_og[u]
        role = 'candidate' if u in strict else ('background' if role == 'non_catalogue_background' and u not in potential and og not in blocked else None)
        key = (species, og, role)
        if key not in needed:
            continue
        labels = sorted(set(cog) & set('ABCDEFGHIJKLMNOPQRSTUVWXYZ') - set('RS'))
        check('missing_COG_not_misread_' + u, cog not in {'NA', 'nan', 'None'})
        counts[key] += 1
        if labels:
            sums[key] += np.array([c in labels for c in CATS], float) / len(labels)
        if role == 'candidate':
            family_units[(species, og)].append(u)
    profiles = {k: sums[k] / counts[k] for k in needed if counts[k]}
    check('all_profiles_available', len(profiles) == len(needed))
    metrics = {(r['species'], r['OG'], r['role']): r for r in rows(F / 'function_robustness/family_annotation_metrics.tsv.gz') if (r['species'], r['OG'], r['role']) in needed}
    check('all_metrics_available', len(metrics) == len(needed))
    for k in needed:
        check('profile_member_count_' + '|'.join(k), counts[k] == int(metrics[k]['units']))
    objects = []
    for (s, o), old in sorted(original.items()):
        controls = pairs[(s, o, 'length_fine_annotation')]
        check('control_number_' + s + o, 3 <= len(controls) <= 5 and len(set(controls)) == len(controls))
        check('background_exclusion_' + s + o, not (set(controls) & blocked))
        p = profiles[(s, o, 'candidate')]
        q = np.stack([profiles[(s, b, 'background')] for b in controls])
        check('family_profile_reproduction_' + s + o, np.allclose(p, old['P'], atol=1e-12, rtol=0) and np.allclose(q.mean(0), old['B'], atol=1e-12, rtol=0))
        check('strict_units_' + s + o, len(family_units[(s, o)]) == old['units'])
        alternate = pairs.get((s, o, 'length_only'), [])
        objects.append(dict(species=s, OG=o, P=p, Q=q, B=q.mean(0), D=p-q.mean(0), controls=controls,
                            units=old['units'], unit_ids=family_units[(s, o)], alternate=alternate,
                            alternate_B=np.mean([profiles[(s, b, 'background')] for b in alternate], axis=0) if alternate else None))
    for r in rows(F / 'ecology/species_function_atlas_data.tsv'):
        if r['species'] in SP:
            obs = [o for o in objects if o['species'] == r['species']]
            j = CATS.index(r['category'])
            check('species_effect_reproduction_' + r['species'] + r['category'], abs(np.mean([o['D'][j] for o in obs])-float(r['effect'])) < 1e-12)
    return objects, strict, blocked, profiles, metrics


def taxonomy():
    parent, phyla = {}, set()
    INPUTS.update([TAX / 'nodes.dmp', TAX / 'names.dmp'])
    with (TAX / 'nodes.dmp').open() as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')[:3]]
            t, p = int(parts[0]), int(parts[1])
            parent[t] = p
            if parts[2] == 'phylum':
                phyla.add(t)
    names, aliases = {}, defaultdict(set)
    with (TAX / 'names.dmp').open() as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')[:4]]
            t = int(parts[0])
            if t in phyla:
                aliases[parts[1]].add(t)
                if parts[3] == 'scientific name':
                    names[t] = parts[1]
    def lineage(t):
        seen = set()
        while t in parent and t not in seen:
            seen.add(t)
            t = parent[t]
        return seen
    bacterial_phyla = {p for p in phyla if 2 in lineage(p)}
    cache = {}
    def resolve(t, saved, allow_saved):
        if t not in cache:
            ancestors = lineage(t)
            ps = ancestors & bacterial_phyla
            cache[t] = next(iter(ps)) if 2 in ancestors and len(ps) == 1 else None
        if cache[t] is not None:
            return cache[t], 'taxid_in_20230804_snapshot'
        if allow_saved:
            ps = set()
            for token in saved.split(';'):
                ps.update(aliases.get(token.strip(), set()) & bacterial_phyla)
            if len(ps) == 1 and 'Bacteria' in saved.split(';'):
                return next(iter(ps)), 'unchanged_reference_saved_lineage_phylum'
        return None, 'no_unambiguous_bacterial_phylum'
    return resolve, names


def add_sources(objects, strict):
    catalogue = {r['hgt_id']: r for r in rows(F / 'inputs/final_catalogue_8989.tsv')}
    check('final_catalogue_8989', len(catalogue) == 8989)
    needed = {r for ob in objects for u in ob['unit_ids'] for r in splitids(strict[u]['record_ids'])}
    check('all_selected_records_in_final_catalogue', needed <= set(catalogue))
    resolve, names = taxonomy()
    sources = {}
    for rid in sorted(needed):
        r = catalogue[rid]
        domain = r['display_best_distant_domain']
        tx = integer(r['display_best_distant_taxid'])
        original_tx = integer(r['Target_taxid'])
        review = r['review_source_annotation_after_model_review']
        phylum, method = None, ''
        if review in {'source_unresolved_after_expanded_sampling', 'quality_flag_prevents_source_revision', 'no_evaluable_rebuilt_tree'}:
            status = 'existing_source_review_unresolved'
        elif domain != 'Bacteria':
            status = 'non_bacterial_or_unresolved_reference'
        else:
            phylum, method = resolve(tx, r['Target_lineage'], tx is not None and tx == original_tx)
            status = 'resolved_bacterial_reference' if phylum else 'unresolved_phylum'
        sources[rid] = dict(record_id=rid, species=r['Specie'], status=status, phylum_taxid=phylum,
            phylum=names.get(phylum, ''), assignment_method=method, display_domain=domain,
            display_taxid=r['display_best_distant_taxid'], original_taxid=r['Target_taxid'],
            display_basis=r['display_best_distant_basis'], existing_source_review=review,
            original_manuscript_category=r['manuscript_HGT_category'], original_evidence_tier=r['evidence_tier'],
            viral_homology=truth(r['viral_annotation_viral_homology_detected']), original_long_branch=r['long_branch'])
    write('record_source_correspondence.tsv', list(sources.values()))
    for ob in objects:
        records = sorted({r for u in ob['unit_ids'] for r in splitids(strict[u]['record_ids'])})
        ss = [sources[r] for r in records]
        ps = {s['phylum_taxid'] for s in ss if s['status'] == 'resolved_bacterial_reference'}
        if not records:
            state = 'no_record_correspondence'
        elif any(s['status'] != 'resolved_bacterial_reference' for s in ss):
            state = 'contains_unresolved_or_non_bacterial_reference'
        elif len(ps) != 1:
            state = 'multiple_bacterial_phyla'
        else:
            state = 'single_consistent_bacterial_phylum'
        ob['source_status'] = state
        ob['source'] = names[next(iter(ps))] if state == 'single_consistent_bacterial_phylum' else None
        ob['records'] = records
        ob['record_source_states'] = ';'.join(sorted({s['status'] for s in ss}))
        ob['all_phyla'] = ';'.join(sorted(names[p] for p in ps))
        ob['viral_homology'] = any(s['viral_homology'] for s in ss)
        ob['tree_classes'] = ';'.join(sorted({s['original_manuscript_category'] for s in ss}))
        ob['evidence_tiers'] = ';'.join(sorted({s['original_evidence_tier'] for s in ss}))
        ob['legacy_no_long_branch_all'] = all(str(s['original_long_branch']).upper() == 'FALSE' for s in ss)
    return catalogue, sources


def weights(obs, stage, target, multiplier=None):
    if not obs:
        return None
    m = np.array([multiplier[o['OG']] for o in obs]) if multiplier is not None else np.ones(len(obs))
    if stage != 'standardized':
        return m/m.sum() if m.sum() else None
    w = np.zeros(len(obs))
    for source, q in target.items():
        ix = np.array([o['source'] == source for o in obs])
        den = m[ix].sum()
        if not den:
            return None
        w[ix] = q*m[ix]/den
    return w


def estimate(obs, stage, target, multiplier=None, ranks=None, alternate=False):
    w = weights(obs, stage, target, multiplier)
    if w is None:
        return None
    p = w @ np.stack([o['P'] for o in obs])
    if alternate:
        if any(o['alternate_B'] is None for o in obs):
            return None
        bg = np.stack([o['alternate_B'] for o in obs])
    elif ranks is not None:
        bg = np.stack([o['Q'][np.argmin([ranks[b] for b in o['controls']])] for o in obs])
    else:
        bg = np.stack([o['B'] for o in obs])
    b = w @ bg
    return p, b, p-b


def common_target(objects, minimum):
    counts = {s: Counter(o['source'] for o in objects if o['species'] == s and o['source']) for s in SP}
    common = sorted(p for p in set(counts[SP[0]]) & set(counts[SP[1]]) if all(counts[s][p] >= minimum for s in SP))
    den = {s: sum(counts[s][p] for p in common) for s in SP}
    target = {p: sum(counts[s][p]/den[s] for s in SP)/2 for p in common} if common else {}
    return counts, common, target


def quantile_rows(storage, kind):
    out = []
    for (stage, comparison), values in storage.items():
        vv = np.stack(values)
        for j, c in enumerate(CATS):
            out.append(dict(kind=kind, stage=stage, comparison=comparison, category=c, draws=len(vv),
                q025_pp=100*np.quantile(vv[:, j], .025), q975_pp=100*np.quantile(vv[:, j], .975),
                median_pp=100*np.median(vv[:, j]), formal_pvalue='NA'))
    return out


def main():
    INPUTS.update([Path(__file__), B/'PROTOCOL_zh.md'])
    objects, strict, blocked, profiles, metrics = read_profiles()
    print('Original 321/263-family profiles reproduced.', flush=True)
    catalogue, record_sources = add_sources(objects, strict)
    counts, common, target = common_target(objects, 5)
    check('main_common_sources_exist', len(common) >= 1)
    panels = {stage: {s: [o for o in objects if o['species'] == s and (stage == 'original' or o['source'] in common)] for s in SP} for stage in STAGES}
    write('source_counts.tsv', [dict(species=s, source=p, families=counts[s][p], primary_common=p in common,
        primary_target_fraction=target.get(p, np.nan)) for s in SP for p in sorted(set(counts[SP[0]]) | set(counts[SP[1]]))])
    write('source_status_counts.tsv', [dict(species=s, source_status=st, families=n) for s in SP for st,n in Counter(o['source_status'] for o in objects if o['species']==s).items()])
    coverage, estimates, contrasts, weightrows, unitrows, pairrows, annotation = [], [], [], [], [], [], []
    point = {}
    for stage in STAGES:
        point[stage] = {}
        for s in SP:
            obs = panels[stage][s]
            w = weights(obs, stage, target)
            p,b,d = estimate(obs, stage, target)
            point[stage][s] = d
            bw = Counter()
            for ob, wt in zip(obs, w):
                for bg in ob['controls']:
                    bw[bg] += wt/len(ob['controls'])
                weightrows.append(dict(stage=stage,species=s,OG=ob['OG'],source=ob['source'],weight=wt,
                    weight_relative_to_common_equal=wt*len(obs)))
            check('normalized_weights_'+stage+s, abs(w.sum()-1)<1e-12 and abs(sum(bw.values())-1)<1e-12)
            coverage.append(dict(stage=stage,species=s,families=len(obs),retained_fraction=len(obs)/sum(o['species']==s for o in objects),
                strict_loci=sum(o['units'] for o in obs),record_ids=len({r for o in obs for r in o['records']}),
                source_categories=len({o['source'] for o in obs if o['source']}),unique_background_OGs=len(bw),
                candidate_weight_ESS=1/sum(w*w),background_weight_ESS=1/sum(x*x for x in bw.values()),
                max_candidate_weight=max(w),max_background_weight=max(bw.values()),viral_homology_families=sum(o['viral_homology'] for o in obs)))
            for j,c in enumerate(CATS):
                estimates.append(dict(stage=stage,species=s,category=c,families=len(obs),candidate_share=p[j],background_share=b[j],effect=d[j],effect_pp=100*d[j]))
            for field in ['eggnog_fraction','informative_cog_fraction','ko_fraction','pfam_fraction','stronger_match_fraction']:
                cv=np.array([float(metrics[(s,o['OG'],'candidate')][field]) for o in obs])
                bv=np.array([np.mean([float(metrics[(s,bg,'background')][field]) for bg in o['controls']]) for o in obs])
                annotation.append(dict(stage=stage,species=s,metric=field,candidate=w@cv,background=w@bv,difference=w@(cv-bv)))
        for j,c in enumerate(CATS):
            contrasts.append(dict(stage=stage,category=c,contrast='T.foetus_minus_T.vaginalis',effect_difference_pp=100*(point[stage][SP[1]][j]-point[stage][SP[0]][j]),formal_pvalue='NA'))
    for ob in objects:
        for j,c in enumerate(CATS):
            unitrows.append(dict(species=ob['species'],OG=ob['OG'],category=c,candidate_share=ob['P'][j],background_share=ob['B'][j],difference=ob['D'][j],
                candidate_loci=ob['units'],source_status=ob['source_status'],source=ob['source'],all_phyla=ob['all_phyla'],record_source_states=ob['record_source_states'],
                primary_common=ob['source'] in common,viral_homology=ob['viral_homology'],legacy_no_long_branch_all=ob['legacy_no_long_branch_all'],
                original_manuscript_categories=ob['tree_classes'],original_evidence_tiers=ob['evidence_tiers'],record_ids=';'.join(ob['records']),unit_ids=';'.join(ob['unit_ids'])))
        for k,bg in enumerate(ob['controls']):
            for j,c in enumerate(CATS):
                pairrows.append(dict(species=ob['species'],candidate_OG=ob['OG'],background_OG=bg,category=c,background_share=ob['Q'][k,j]))
    write('family_function_source.tsv',unitrows)
    write('original_matched_backgrounds.tsv',pairrows)
    write('coverage.tsv',coverage)
    write('family_weights.tsv',weightrows)
    write('species_effects.tsv',estimates)
    write('between_species_contrasts.tsv',contrasts)
    write('annotation_balance.tsv',annotation)
    print('Common source categories:', len(common), 'families:', [len(panels['common'][s]) for s in SP], flush=True)

    rng=np.random.default_rng(SEED)
    all_ogs=sorted({o['OG'] for o in objects})
    samples=defaultdict(list)
    rejects=Counter()
    for rep in range(1000):
        mult=dict(zip(all_ogs,rng.poisson(1,len(all_ogs))))
        for stage in STAGES:
            ee={s:estimate(panels[stage][s],stage,target,multiplier=mult) for s in SP}
            if any(v is None for v in ee.values()):
                rejects[stage]+=1
                continue
            for s in SP:
                samples[(stage,s)].append(ee[s][2])
            samples[(stage,'T.foetus_minus_T.vaginalis')].append(ee[SP[1]][2]-ee[SP[0]][2])
    intervals=quantile_rows(samples,'shared_OG_weight_sensitivity')
    all_bg=sorted({b for o in objects for b in o['controls']})
    draws=defaultdict(list)
    for rep in range(500):
        ranks=dict(zip(all_bg,rng.random(len(all_bg))))
        for stage in STAGES:
            ee={s:estimate(panels[stage][s],stage,target,ranks=ranks) for s in SP}
            for s in SP:
                draws[(stage,s)].append(ee[s][2])
            draws[(stage,'T.foetus_minus_T.vaginalis')].append(ee[SP[1]][2]-ee[SP[0]][2])
    intervals += quantile_rows(draws,'fixed_matched_background_choice_sensitivity')
    write('descriptive_sensitivity_ranges.tsv',intervals)

    totals=Counter()
    for o in objects:
        totals[o['OG']]+=o['units']
    largest=set(sorted(totals,key=lambda o:(-totals[o],o))[:5])
    shared=set(o['OG'] for o in objects if o['species']==SP[0]) & set(o['OG'] for o in objects if o['species']==SP[1])
    write('fixed_exclusion_sets.tsv',[dict(OG=o,total_selected_candidate_loci=totals[o],shared_candidate_OG=o in shared,top5_locus_count=o in largest) for o in sorted(shared|largest)])
    sensitivity=[]
    for name,excluded in [('remove_shared_candidate_OGs',shared),('remove_top5_locus_count_OGs',largest)]:
        for stage in STAGES:
            oo={s:[o for o in panels[stage][s] if o['OG'] not in excluded] for s in SP}
            ee={s:estimate(oo[s],stage,target) for s in SP}
            for j,c in enumerate(CATS):
                sensitivity.append(dict(scenario=name,stage=stage,category=c,families_Tv=len(oo[SP[0]]),families_Tf=len(oo[SP[1]]),
                    effect_Tv_pp=100*ee[SP[0]][2][j] if ee[SP[0]] is not None else np.nan,
                    effect_Tf_pp=100*ee[SP[1]][2][j] if ee[SP[1]] is not None else np.nan,
                    contrast_pp=100*(ee[SP[1]][2][j]-ee[SP[0]][2][j]) if all(v is not None for v in ee.values()) else np.nan,
                    estimable=all(v is not None for v in ee.values()),target_weights='main_fixed'))
    _,allcommon,alltarget=common_target(objects,1)
    for name,alternate,cset,tt in [('all_shared_phyla_min1',False,allcommon,alltarget),('existing_length_only_controls_same_candidates',True,common,target)]:
        for stage in STAGES:
            oo={s:[o for o in objects if o['species']==s and (stage=='original' or o['source'] in cset)] for s in SP}
            ee={s:estimate(oo[s],stage,tt,alternate=alternate) for s in SP}
            for j,c in enumerate(CATS):
                sensitivity.append(dict(scenario=name,stage=stage,category=c,families_Tv=len(oo[SP[0]]),families_Tf=len(oo[SP[1]]),
                    effect_Tv_pp=100*ee[SP[0]][2][j] if ee[SP[0]] is not None else np.nan,
                    effect_Tf_pp=100*ee[SP[1]][2][j] if ee[SP[1]] is not None else np.nan,
                    contrast_pp=100*(ee[SP[1]][2][j]-ee[SP[0]][2][j]) if all(v is not None for v in ee.values()) else np.nan,
                    estimable=all(v is not None for v in ee.values()),target_weights='all_shared_equal_species' if name.startswith('all_shared') else 'main_fixed'))
    write('prespecified_sensitivities.tsv',sensitivity)
    deletion=[]
    for stage in STAGES:
        active=sorted({o['OG'] for s in SP for o in panels[stage][s]})
        vv=[]
        for omitted in active:
            oo={s:[o for o in panels[stage][s] if o['OG'] != omitted] for s in SP}
            ee={s:estimate(oo[s],stage,target) for s in SP}
            if all(v is not None for v in ee.values()):
                vv.append(ee[SP[1]][2]-ee[SP[0]][2])
        a=np.stack(vv)
        for j,c in enumerate(CATS):
            deletion.append(dict(stage=stage,category=c,valid_deletions=len(a),minimum_contrast_pp=100*a[:,j].min(),maximum_contrast_pp=100*a[:,j].max(),
                any_sign_change=bool(np.any(np.sign(a[:,j]) != np.sign(point[stage][SP[1]][j]-point[stage][SP[0]][j])))))
    write('leave_one_OG_sensitivity.tsv',deletion)
    support=[]
    for stage in STAGES:
        for s in SP:
            obs=panels[stage][s]
            support.append(dict(stage=stage,species=s,families=len(obs),all_legacy_no_long_branch=sum(o['legacy_no_long_branch_all'] for o in obs),
                any_donor_leaning_or_placed=sum('Donor-' in o['tree_classes'] for o in obs),
                lineage_retained_only=sum(o['tree_classes']=='Lineage-retained ancestral HGT candidates' for o in obs),
                viral_homology=sum(o['viral_homology'] for o in obs),source_direction_newly_confirmed=0))
    write('preserved_evidence_coverage.tsv',support)
    fig_effects(estimates,contrasts,intervals)
    fig_sources(counts,common,target,coverage)
    write('validation_checks.tsv',CHECKS)
    write('input_manifest.tsv',[dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(INPUTS)])
    status=dict(status='completed_descriptive_composition_analysis',species=SP,original_families=[321,263],common_sources=common,
        common_families=[len(panels['common'][s]) for s in SP],all_shared_sources=len(allcommon),shared_candidate_OGs=len(shared),
        requested_OG_draws=1000,rejected_zero_stratum_draws=dict(rejects),background_draws=500,
        checks_passed=len(CHECKS),formal_pvalues_computed=False,source_direction_reassigned=False,original_catalogue_modified=False,
        safety_skipped=[],new_sequence_analysis=False)
    (B/'completion.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(status,ensure_ascii=False),flush=True)


def fig_effects(estimates,contrasts,intervals):
    df=pd.DataFrame(estimates);dd=pd.DataFrame(contrasts)
    it=pd.DataFrame(intervals)
    fig,axes=plt.subplots(2,4,figsize=(15,7),sharex=True)
    colors=['#2b6f9c','#bf6434']
    titles=['E: amino acid','F: nucleotide','G: carbohydrate','H: coenzyme']
    for j,c in enumerate(CATS):
        ax=axes[0,j]
        for s,short,color,offset in zip(SP,SHORT,colors,[-.04,.04]):
            vals=[float(df[(df.species==s)&(df.category==c)&(df.stage==st)].effect_pp.iloc[0]) for st in STAGES]
            ax.plot(np.arange(3)+offset,vals,'o-',color=color,label=short,lw=1.5)
        ax.axhline(0,color='#aaaaaa',lw=.8);ax.set_title(titles[j]);ax.grid(axis='y',alpha=.2)
        bot=axes[1,j]
        vals=[float(dd[(dd.category==c)&(dd.stage==st)].effect_difference_pp.iloc[0]) for st in STAGES]
        lo=[];hi=[]
        for st in STAGES:
            z=it[(it.kind=='shared_OG_weight_sensitivity')&(it.comparison=='T.foetus_minus_T.vaginalis')&(it.category==c)&(it.stage==st)].iloc[0]
            lo.append(z.q025_pp);hi.append(z.q975_pp)
        bot.vlines(np.arange(3),lo,hi,color='#777777',lw=2)
        bot.plot(vals,'o-',color='#292929');bot.axhline(0,color='#aaaaaa',lw=.8);bot.grid(axis='y',alpha=.2)
        bot.set_xticks(range(3),['Original','Common','Standardized'],rotation=20)
    axes[0,0].set_ylabel('Candidate minus background (pp)')
    axes[1,0].set_ylabel('T. foetus minus T. vaginalis (pp)')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Two-species source-conditioned functional composition',fontsize=15)
    fig.text(.5,.012,'Vertical ranges: 2.5–97.5% shared-OG weight sensitivity; not phylogenetic confidence intervals or ecological tests.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.055,1,.95])
    for ext in ['png','pdf','svg']:
        fig.savefig(B/('source_function_effects.'+ext),dpi=180)
    plt.close(fig)


def fig_sources(counts,common,target,coverage):
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    x=np.arange(len(common));width=.25
    for s,offset,label,color in [(SP[0],-width,SHORT[0],'#2b6f9c'),(SP[1],0,SHORT[1],'#bf6434')]:
        n=sum(counts[s][p] for p in common)
        axes[0].bar(x+offset,[counts[s][p]/n for p in common],width,label=label,color=color)
    axes[0].bar(x+width,[target[p] for p in common],width,label='Equal-species target',color='#909090')
    axes[0].set_xticks(x,common,rotation=40,ha='right');axes[0].set_ylabel('Fraction within common range');axes[0].legend(fontsize=8)
    axes[0].set_title('Source affinity, not inferred transfer direction')
    for i,s in enumerate(SP):
        base=next(r for r in coverage if r['stage']=='original' and r['species']==s)
        sub=next(r for r in coverage if r['stage']=='common' and r['species']==s)
        axes[1].bar(i,base['families'],color='#dddddd')
        axes[1].bar(i,sub['families'],color=['#2b6f9c','#bf6434'][i])
        axes[1].text(i,base['families']+4,str(base['families']),ha='center')
        axes[1].text(i,sub['families']/2,f"{sub['families']}\n{sub['retained_fraction']:.1%}",ha='center',color='white')
    axes[1].set_xticks([0,1],SHORT);axes[1].set_ylabel('Candidate species × OG units');axes[1].set_title('Common-range coverage (coloured)')
    axes[1].set_ylim(0,360)
    fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(B/('source_coverage.'+ext),dpi=180)
    plt.close(fig)


if __name__=='__main__':
    main()
