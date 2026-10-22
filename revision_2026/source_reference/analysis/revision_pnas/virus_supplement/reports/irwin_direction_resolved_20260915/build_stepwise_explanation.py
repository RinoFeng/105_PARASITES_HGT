from pathlib import Path
import pandas as pd,json
R=Path(__file__).resolve().parent
T=pd.read_csv(R/'combined_1013_traced_historical_matches.tsv.gz',sep='\t',keep_default_na=False)
V=pd.read_csv(R/'final_event_historical_protein_filters.tsv.gz',sep='\t',keep_default_na=False)
stages=[('完整序列可比集合',None),('原始分布／相似性／熵初筛','in_saved_entropy_pass'),('BUSCO保守同源过滤','in_post_busco'),('模式动物保守同源过滤','in_post_conserved'),('基因组存在性过滤','in_post_genome_presence'),('EGT／共生来源过滤','in_post_egt_clean'),('早期映射／contig验证','in_post_contig_validation'),('种内代表选择及整合','in_integrated'),('修订contig分级','in_tree_targets'),('取得可评估树','tree_primary_class'),('原树分类保留名单','in_final_catalog')]
rows=[];excluded=[]
for direction in ['VtoE','unkn']:
 ids=set(V.loc[(V.direction==direction)&V.archive_exact_local,'irwin_table2_id']);a=T[T.irwin_table2_id.isin(ids)].set_index('irwin_table2_id');prev=ids
 for i,(name,col) in enumerate(stages):
  keep=ids if col is None else set(a.loc[a[col]!=''].index) if col=='tree_primary_class' else set(a.loc[a[col].astype(str).str.lower()=='true'].index)
  assert keep.issubset(prev),(direction,name,keep-prev)
  removed=prev-keep;rows.append(dict(direction=direction,stage_order=i,stage=name,entering=len(prev),removed=len(removed),remaining=len(keep)))
  for sid in sorted(removed):
   q=a.loc[sid];excluded.append(dict(direction=direction,irwin_table2_id=sid,accession=q.accession,genome=q.genome,stage=name,terminal_reason=q.historical_terminal_stage,initial_first_failure=q.old_replay_first_failure,all_initial_failures=q.old_replay_all_failures,tree_class=q.tree_primary_class,contig_tier=q.contig_context_tier,cluster_representative=q.intraspecies_cluster_representative,new_terminal_stage=q.new_terminal_stage,new_followup_all_failures=q.new_followup_all_failures))
  prev=keep
f=pd.DataFrame(rows);d=pd.DataFrame(excluded);f.to_csv(R/'verified_stepwise_original_funnel.tsv',sep='\t',index=False);d.to_csv(R/'stepwise_removed_protein_list.tsv',sep='\t',index=False)
follow=V[V.direction.isin(['VtoE','unkn'])&(V.new_terminal_stage=='followup_locus_not_retained')][['direction','irwin_table2_id','archive_local_accessions','new_followup_all_failures']].drop_duplicates();assert len(follow)==11;follow.to_csv(R/'11_viral_followup_rejections.tsv',sep='\t',index=False)
with pd.ExcelWriter(R/'未列入事件解释与逐步过滤.xlsx') as out:
 f.to_excel(out,sheet_name='原流程逐步漏斗',index=False);d.to_excel(out,sheet_name='逐条退出记录',index=False);follow.to_excel(out,sheet_name='病毒补筛后续11条',index=False);pd.read_csv(R/'15_unlisted_overlap_family_tree_audit.tsv',sep='\t').to_excel(out,sheet_name='原15条事件成员审计',index=False)
assert sum(f.loc[(f.direction=='VtoE'),'removed'])==537 and sum(f.loc[(f.direction=='unkn'),'removed'])==364
(R/'stepwise_validation.json').write_text(json.dumps({'all_stage_sets_nested':True,'VtoE_start':544,'VtoE_original_retained':7,'unkn_start':385,'unkn_original_retained':21,'mutually_exclusive_removed_records':len(d),'viral_followup_unretained':len(follow)},indent=2)+'\n')
print(f.to_string(index=False))
