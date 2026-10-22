from pathlib import Path
import csv,json,collections,hashlib
D=Path(__file__).resolve().parents[1];B=D.parent;A=B/'harmonized_final_20260915';T=D/'tables'
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(name,rr):
 fields=list(dict.fromkeys(k for r in rr for k in r))
 with (T/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def yes(s):return str(s).lower()=='true'
def domain(s):
 for k in ['Viruses','Bacteria','Archaea','Eukaryota']:
  if k in s.split(';'):return k
 return 'unresolved'
old=rows(A/'inputs/historical8985_frozen.tsv');new=rows(A/'inputs/discovery120_frozen.tsv');audit={x['locus_id']:x for x in rows(A/'reports/discovery120_original_rule_audit.tsv')};ev={x['hgt_id']:x for x in rows(B/'reports/existing_catalog_viral_evidence_with_original_function.tsv')};rank={x['hgt_id']:x for x in rows(B/'reports/merged_viral_ranking/catalogue_8985_merged_ranking_summary.tsv')};rev={x['hgt_id']:x for x in rows(B/'reports/donor_reassessment/original74_donor_revision_proposals.tsv')};frac={x['hgt_id']:x for x in rows(T/'old8985_viral_out_fractions.tsv')};context={x['locus_id']:x for x in rows(B/'reports/novel_loci_retained_with_bilateral_flanks.tsv')}
selected=[x for x in new if x['viral_vs_cellular_out'] in ['viral_higher','no_cellular_reference']];selected_ids={x['locus_id'] for x in selected}
ir=rows(B/'reports/irwin_direction_resolved_20260915/all54_S2_catalogue_sequence_overlaps.tsv');bykey=collections.defaultdict(list)
for x in ir:
 for g in x['archive_local_genomes'].split(';'):
  for a in x['archive_local_accessions'].split(';'):
   if g and a:bykey[g,a].append(x)
def literature(g,a):
 z=bykey.get((g,a),[]);primary=[x for x in z if float(x['contamination_score'])>=2];explicit=[x for x in primary if x['final_event_ids']]
 return dict(Irwin_full_sequence_exact_ids=';'.join(sorted({x['irwin_table2_id'] for x in z})),Irwin_primary_score_ge2_ids=';'.join(sorted({x['irwin_table2_id'] for x in primary})),Irwin_final_event_ids=';'.join(sorted({e for x in explicit for e in x['final_event_ids'].split(';') if e})),Irwin_final_directions=';'.join(sorted({x['direction'] for x in explicit if x['direction']})),Irwin_unlisted_final_member_ids=';'.join(sorted({x['irwin_table2_id'] for x in primary if not x['final_event_ids']})),Irwin_low_score_ids=';'.join(sorted({x['irwin_table2_id'] for x in z if float(x['contamination_score'])<2})),Irwin_matching_method='archived_full_amino_acid_sequence_equality;accession_aliases_allowed')
# Resolve domain only for already reported best taxids, preserving all original donor fields.
needed={x['uniform_best_cellular_out_taxid'] for x in rev.values()};taxdomains={}
root=A.parents[3]
for p in [root/'align/protein_align/genome_divide/target_nt_taxid_lineage.txt',root/'align/nucl_align/genome_divided/target_nt_taxid_lineage.txt',root/'align/target_nt_taxid_lineage2.txt']:
 with p.open() as f:
  for l in f:
   z=l.rstrip('\n').split('\t')
   if z[0] in needed and len(z)>1:taxdomains[z[0]]=domain(z[1])
primary=[];plot=[];decisions=[]
for x in old:
 h=x['hgt_id'];r=rank[h];v=rev.get(h);out=dict(x);out.update(cohort_role='primary_historical8985',counting_unit='historical_catalogue_record',primary_plot_include=True,final_addition=False)
 out.update({'viral_annotation_'+k:val for k,val in ev.get(h,{}).items() if k not in ['hgt_id','species','genome','accession','data_type']})
 out.update({'saved_rank_'+k:val for k,val in r.items() if k not in ['hgt_id','species','genome','accession','data_type']})
 out.update({'review_'+k:val for k,val in (v or {}).items() if k not in ['hgt_id','species','accession','data_type']})
 out.update({k:val for k,val in frac[h].items() if k not in ['hgt_id','species','data_type']});out.update(literature(x['Genome'],x['Accession']))
 status=v['uniform_viral_vs_cellular_out'] if v else 'not_uniformly_researched';best=x['Target_accesson'];tid=x['Target_taxid'];name=x['Target_lineage'].split(';')[-1];dom=domain(x['Target_lineage']);basis='original_reported_distant_reference_not_uniformly_reassessed'
 if v and status=='viral_higher':best=v['uniform_best_viral'];tid=v['uniform_best_viral_taxid'];name=v['uniform_best_viral_taxon_name'];dom='Viruses';basis='highest_classified_distant_bitscore_in_uniform_search'
 elif v and status in ['cellular_higher','no_viral_hit'] and v['uniform_best_cellular_out']:
  best=v['uniform_best_cellular_out'];tid=v['uniform_best_cellular_out_taxid'];name=v['uniform_best_cellular_out_taxon_name'];dom=taxdomains.get(tid,'cellular_domain_unresolved');basis='highest_classified_distant_bitscore_in_uniform_search'
 elif v:basis='uniform_search_has_no_cellular_out_reference;historical_reference_preserved'
 out.update(display_best_distant_accession=best,display_best_distant_taxid=tid,display_best_distant_taxon=name,display_best_distant_domain=dom,display_best_distant_basis=basis,uniform_viral_rank_status=status,uniform_search_reviewed=bool(v),uniform_unknown_competes=v['unknown_competes'] if v else 'not_assessed',historical_donor_direction_changed=False,documented_new_tree_conflict=bool(v and v['source_annotation_after_model_review']=='tree_class_outside_retained_whitelist'))
 primary.append(out)
 plot.append(dict(record_id=h,cohort='primary_historical8985',primary_plot_include=True,species=x['Specie'],genome=x['Genome'],accession=x['Accession'],data_type=x['Data_type'],historical_donor_domain=domain(x['Target_lineage']),display_best_distant_domain=dom,display_best_distant_taxon=name,display_best_distant_basis=basis,uniform_viral_rank_status=status,current_viral_panel_positive=r['current_viral_hit_present'],historical_tree_category=x['manuscript_HGT_category'],documented_new_tree_conflict=out['documented_new_tree_conflict'],contig_context_tier=x['contig_context_tier'],COG_category=x['COG_category'],KEGG_ko=x['KEGG_ko'],Description=x['Description'],GC_Content_historical=x['GC_Content'],length_historical=x['Length'],length_historical_data_type=x['Data_type'],**literature(x['Genome'],x['Accession'])))
novel=[]
for x in new:
 a=audit[x['locus_id']];c=context[x['locus_id']];keep=x['locus_id'] in selected_ids
 d=dict(x);d.update({'audit_'+k:v for k,v in a.items() if k not in x});d.update({'context_'+k:v for k,v in c.items() if k not in x});d.update(literature(x['genome'],x['accession']))
 conflicts=[v for v in a['blocking_or_dedup_flags'].split(';') if v and v!='classification_only_initial_not_passed'];pending=['final_unified_eligible_pool_dedup','original_tree_sampling_equivalence']
 if a['classification_only_initial_pass']!='True':pending.append('database_addition_branch_initial_rules_not_equivalent')
 if x['workflow']=='nucleotide':pending.extend(['identical_region_initial_metrics','original_nucleotide_followup_chain'])
 d.update(cohort_role='viral_priority_supplement_pending' if keep else 'outside_viral_priority_reporting_scope',primary_plot_include=False,viral_priority_selected=keep,selection_rule='uniform_viral_higher_or_no_qualifying_cellular_out',selection_design='exploratory_after_existing_results',final_addition_status='not_certified',documented_followup_conflicts=';'.join(conflicts),remaining_checks=';'.join(pending),not_cellular_out_does_not_mean_no_self_homologs=True)
 decisions.append(dict(locus_id=x['locus_id'],locus_cluster=x['locus_cluster'],accession=x['accession'],species=x['species'],viral_vs_cellular_out=x['viral_vs_cellular_out'],selected=keep,reason='viral_priority_scope' if keep else 'cellular_out_scores_higher;not_a_claim_of_no_HGT',documented_followup_conflicts=';'.join(conflicts),pending_checks=';'.join(pending),final_addition_status='not_certified'))
 if keep:
  novel.append(d);plot.append(dict(record_id=x['locus_id'],cohort='viral_priority_supplement_pending',primary_plot_include=False,species=x['species'],genome=x['genome'],accession=x['accession'],data_type=x['workflow'],locus_cluster=x['locus_cluster'],display_best_distant_domain='Viruses',display_best_distant_taxon=x['best_viral_taxon_label'],display_best_distant_basis='viral_higher_than_cellular_out' if x['viral_vs_cellular_out']=='viral_higher' else 'viral_detected_without_qualifying_cellular_out_reference',uniform_viral_rank_status=x['viral_vs_cellular_out'],current_viral_panel_positive=True,historical_tree_category=c['manuscript_tree_category'],documented_followup_conflicts=';'.join(conflicts),tree_flags=c['tree_flags'],region_GC_percent=x['region_GC_percent'],query_aa_length=x['query_aa_length_actual'],genomic_span_bp=x['genomic_span_bp'],minimal_HGT_age_MRCA_Ma=x['raw_MRCA_Ma'],minimal_age_interpretation=x['interpretation'],**literature(x['genome'],x['accession'])))
write('primary8985_with_viral_annotation.tsv',primary);write('viral_priority19_all_evidence.tsv',novel);write('all120_scope_decisions.tsv',decisions);write('plotting_records_with_cohort_flags.tsv',plot)
clusters=[]
for k in sorted({x['locus_cluster'] for x in novel}):
 z=[x for x in novel if x['locus_cluster']==k]
 clusters.append(dict(locus_cluster=k,species=z[0]['species'],genome=z[0]['genome'],records=';'.join(x['locus_id'] for x in z),accessions=';'.join(sorted({x['accession'] for x in z})),record_count=len(z),ranking_states=';'.join(sorted({x['viral_vs_cellular_out'] for x in z})),any_documented_conflict=any(bool(x['documented_followup_conflicts']) for x in z),primary_plot_include=False,final_addition_status='not_certified'))
write('viral_priority17_historical_locus_units.tsv',clusters)
comparison=[]
for x in ir:
 key={(g,a) for g in x['archive_local_genomes'].split(';') for a in x['archive_local_accessions'].split(';')}
 oldmatch=[r['hgt_id'] for r in old if (r['Genome'],r['Accession']) in key];newmatch=[r['locus_id'] for r in selected if (r['genome'],r['accession']) in key]
 if oldmatch or newmatch:comparison.append(dict(irwin_protein=x['irwin_table2_id'],contamination_score=x['contamination_score'],primary_literature_comparison=float(x['contamination_score'])>=2,direction=x['direction'] or 'no_explicit_final_event',event_ids=x['final_event_ids'],primary8985_matches=';'.join(oldmatch),viral_priority19_matches=';'.join(newmatch),matching_basis='archived_full_sequence_identity',novel_inclusion_status='supplement_pending_not_final_addition' if newmatch else 'not_applicable'))
write('Irwin_comparison_by_reporting_cohort.tsv',comparison)
summary={}
for label,key in [('primary8985','primary8985_matches'),('viral_priority19','viral_priority19_matches')]:
 z=[x for x in comparison if x[key] and x['primary_literature_comparison']];summary[label]={'proteins_by_direction':dict(collections.Counter(x['direction'] for x in z)),'events_by_direction':{dr:len({e for x in z if x['direction']==dr for e in x['event_ids'].split(';') if e}) for dr in ['VtoE','unkn','EtoV']}}
assert len(primary)==8985 and len(novel)==19 and len(clusters)==17 and len(plot)==9004
status=dict(reporting_main_frozen=True,main_records=8985,original_columns_and_values_preserved=True,viral_priority_records=19,viral_priority_historical_locus_units=17,combined_plotting_rows=9004,combined_plotting_rows_are_not_a_final_HGT_total=True,novel_final_set_certified=False,old_current_viral_panel_positive=sum(yes(rank[x['hgt_id']]['current_viral_hit_present']) for x in old),uniform_reviewed_old_records=len(rev),uniform_viral_higher_old_records=sum(x['uniform_viral_rank_status']=='viral_higher' for x in primary),uniform_viral_higher_among_current1194=sum(x['uniform_viral_rank_status']=='viral_higher' and yes(x['current_viral_panel_positive']) for x in primary),unknown_competition_among_viral_higher=sum(x['uniform_viral_rank_status']=='viral_higher' and yes(x['uniform_unknown_competes']) for x in primary),selected_ranking_counts=dict(collections.Counter(x['viral_vs_cellular_out'] for x in novel)),Irwin=summary,scope_excluded_records=101,no_transfer_vector_inference=True,safety_skipped=[])
(D/'reporting_manifest.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n');print(json.dumps(status,ensure_ascii=False,indent=2))
