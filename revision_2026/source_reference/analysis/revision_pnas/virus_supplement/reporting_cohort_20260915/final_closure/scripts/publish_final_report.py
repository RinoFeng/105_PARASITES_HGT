from pathlib import Path
from collections import defaultdict,Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import csv,json,hashlib
F=Path(__file__).resolve().parents[1];D=F.parent;B=D.parent;S=F/'diamond_review'
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rs):
 if not rs:return
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(dict.fromkeys(k for r in rs for k in r)),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def ids(s):return set(filter(None,(s or '').split(';')))
g=json.load((F/'final_gate_summary.json').open());gates=rows(F/'tables/final_candidate_gate19.tsv');record={r['locus_id']:r for r in rows(D/'tables/viral_priority19_all_evidence.tsv')};metrics={r['locus_id']:r for r in rows(S/'initial_screen_combined19.tsv')};arms={(r['locus_id'],r['arm']):r for r in rows(S/'initial_screen_three_arms.tsv')};accepted=[record[q] for q in g['accepted_locus_ids']];gate={r['locus_id']:r for r in gates};tax=json.load((S/'resolved_taxonomy.json').open());allhits=rows(S/'all_hits_annotated.tsv');old=rows(D/'tables/primary8985_with_viral_annotation.tsv');source_review={r['locus_id']:r for r in rows(F/'tables/viral_source_topology_review.tsv')};categories=json.load((B/'reports/original_tree_category_map.json').open())['categories']
additions=[];features=[];attribution=[]
for r in accepted:
 q=r['locus_id'];m=metrics[q];c=gate[q];sid=m['best_qualified_viral_subject'];hit=next(h for h in allhits if h['locus_id']==q and h['subject']==sid and h['group']=='VIRAL');domain='VIRAL'
 a={k:'' for k in old[0]};a.update(hgt_id=q,Specie=r['species'],Genome=r['genome'],Accession=r['accession'],Data_type='Protein',Length=r['query_aa_length_actual'],GC_Content=r['CDS_GC_percent'],Contig_Name=r['contig'],hgt_start=r['start_1based'],hgt_end=r['end_1based'],primary_class=c.get('tree_primary_class',''),manuscript_HGT_category=categories.get(c.get('tree_primary_class',''),''),treefile=c.get('tree_primary_file',''),Target_accesson=sid,Target_taxid=hit['taxid'],Target_lineage=hit['lineage'],cohort_role='certified_supplementary_candidate',primary_plot_include='True',display_best_distant_accession=sid,display_best_distant_taxid=hit['taxid'],display_best_distant_taxon=hit['lineage'].split(';')[-1],display_best_distant_domain=domain,display_best_distant_basis='uniform_DIAMOND_qualified_distant_homology;potential_source_not_direction_proof',uniform_viral_rank_status=m['viral_rank_qualified'],uniform_search_reviewed='True',uniform_unknown_competes=m['unknown_competes_with_best_distant'])
 for k in old[0]:
  if k.startswith('Irwin_'):a[k]=r.get(k,'')
 a['supplementary_method']='uniform_DIAMOND_with_frozen_external_species_entropy_cutoff';a['supplementary_source_direction']='potential_viral_source_not_proven_historical_direction';a['supplementary_complete_gate_file']=str(F/'tables/final_candidate_gate19.tsv');additions.append(a)
 x=dict(r);x.update({'final_'+k:v for k,v in m.items()});x.update({'final_gate_'+k:v for k,v in c.items()});x.update({'source_tree_'+k:v for k,v in source_review[q].items()});x['prior_viral_KaKs_subject_matches_current_best']=r.get('viral_KaKs_subject','')==sid;features.append(x)
 no=arms[q,'nr_viral_unreported']['initial_distribution_pass']=='True';nr=arms[q,'nr_viral_reported']['initial_distribution_pass']=='True';attribution.append(dict(locus_id=q,accession=r['accession'],old_MMseqs_classification_repair_initial=r['audit_classification_only_initial_pass'],uniform_nr_viral_unreported_pass=no,uniform_nr_viral_reported_pass=nr,uniform_nr_plus_viral_pass=True,contribution='viral_reporting_recovery_in_uniform_DIAMOND' if not no and nr else 'viral_database_panel_increment_in_uniform_DIAMOND' if not nr else 'detectable_without_viral_reporting_in_uniform_DIAMOND'))
write(F/'tables/accepted_new_candidates.tsv',features);write(F/'tables/accepted_candidate_attribution.tsv',attribution)
combined=old+additions;name=('final' if g['final_count_determined_within_declared_scope'] else 'provisional_certified')+'_catalogue_'+str(len(combined))+'.tsv';write(F/'tables'/name,combined)
# Exact original fields must stay unchanged for every historical record.
assert len(old)==8985 and combined[:8985]==old
with (B/'downstream_novel/queries_all.faa').open() as h:aa={t.split()[0]:s for t,s in SimpleFastaParser(h)}
with (B/'downstream_novel/regions_oriented.fna').open() as h:nt={t.split()[0]:s for t,s in SimpleFastaParser(h)}
for filename,data in [('accepted_new.proteins.faa',aa),('accepted_new.regions.fna',nt)]:
 (F/filename).write_text(''.join('>'+r['locus_id']+' '+r['accession']+'\n'+data[r['locus_id']]+'\n' for r in accepted))
# Recompute author event overlap from explicit expanded members, not family membership.
events={r['event_id']:r for r in rows(B/'reports/irwin_direction_resolved_20260915/final_events_with_S1_correspondence.tsv')};relations=rows(B/'reports/irwin_direction_resolved_20260915/final_event_S2_protein_relations.tsv');pairs=defaultdict(set)
for r in relations:
 if r['passes_Irwin_score_cutoff']=='True':
  assert r['event_id'] in events and r['direction']==events[r['event_id']]['direction'];pairs[r['irwin_table2_id']].add((r['direction'],r['event_id']))
overlap=[];proof=[]
for cohort,rr in [('original8985',old),('accepted_new',accepted),('combined',old+accepted)]:
 lit=set().union(*(ids(r.get('Irwin_primary_score_ge2_ids','')) for r in rr)) if rr else set();listed={i for i in lit if pairs[i]}
 for direction in ['VtoE','EtoV','unkn']:
  pp={i for i in lit if any(d==direction for d,e in pairs[i])};ee={e for i in pp for d,e in pairs[i] if d==direction};overlap.append(dict(cohort=cohort,direction=direction,unique_Irwin_proteins=len(pp),unique_Irwin_events=len(ee),contamination_score_cutoff='>=2',match='complete_amino_acid_sequence_equality_not_length_or_accession_only'))
  for i in sorted(pp):
   for d,e in sorted(pairs[i]):
    if d==direction:proof.append(dict(cohort=cohort,irwin_sequence_id=i,event_id=e,direction=d,source=str(B/'reports/irwin_direction_resolved_20260915/final_event_S2_protein_relations.tsv')))
 overlap.append(dict(cohort=cohort,direction='S2_exact_without_explicit_final_event_membership',unique_Irwin_proteins=len(lit-listed),unique_Irwin_events='',contamination_score_cutoff='>=2',match='S2_sequence_match_alone_is_not_author_final_HGT_evidence'))
write(F/'tables/Irwin_final_overlap.tsv',overlap);write(F/'tables/Irwin_final_overlap_proof.tsv',proof)
# Add current event fields without altering any historical field or value.
for a in combined:
 lit=ids(a.get('Irwin_primary_score_ge2_ids',''))
 for direction in ['VtoE','EtoV','unkn']:
  a['final_review_Irwin_'+direction+'_sequence_ids']=';'.join(sorted(i for i in lit if any(d==direction for d,e in pairs[i])))
  a['final_review_Irwin_'+direction+'_event_ids']=';'.join(sorted({e for i in lit for d,e in pairs[i] if d==direction}))
 a['final_review_Irwin_S2_only_ids']=';'.join(sorted(i for i in lit if not pairs[i]))
 a['final_review_Irwin_evidence_source']='explicit_expanded_final_event_members;score_ge2;full_AA_equality'
write(F/'tables'/name,combined)

rankcounts=Counter(r.get('uniform_viral_rank_status','') for r in old if r.get('uniform_search_reviewed')=='True');old_v=[r for r in old if r.get('uniform_search_reviewed')=='True' and r.get('uniform_viral_rank_status') in ['viral_higher','no_cellular_reference']];old_unknown=sum(r.get('uniform_unknown_competes')=='True' for r in old_v)
summary=dict(g,catalogue_file=str(F/'tables'/name),new_best_qualified_distant_viral=len(additions),historical_uniform_review_rank_counts=dict(rankcounts),historical_reviewed_viral_best_or_no_cellular_reference=len(old_v),historical_viral_best_with_competing_unknown=old_unknown,attribution=dict(Counter(r['contribution'] for r in attribution)),Irwin_overlap=overlap,all_old_fields_unchanged=True,remaining_nt_original_chain='not_independently_completed;two_same_locus_alternate_ORFs_not_counted',safety_skipped=[],final_status='completed_with_declared_scope_and_calibration' if g['final_count_determined_within_declared_scope'] else 'unresolved_candidates_remain')
(F/'final_manifest.json').write_text(json.dumps(summary,indent=2))
n=len(accepted);headline=f'本轮固定的 19 条优先记录完成顺序复核后，可新增 **{n} 条候选**；原 **8,985 条**全部保留，合并为 **{8985+n:,} 条记录**。' if g['final_count_determined_within_declared_scope'] else f'目前完成认证 {n} 条；仍有 {len(g["unresolved_locus_ids"])} 条未闭合，尚不能给出最终新增数。'
lines=['# 病毒候选最终纳入复核报告','',headline,'','这是本轮补充候选在下述明确方法下的结果。原目录没有重新执行本轮全部 DIAMOND 筛选，不能称整个合并目录经过同一次统一重筛；也不把这些记录数解释为独立转移事件数。','','## 纳入条目','','| 蛋白 | 物种 | 统一搜索中的增量来源 | 原模型树类别 | 病毒来源含义 |','|---|---|---|---|---|']
for r in accepted:
 a=next(x for x in attribution if x['locus_id']==r['locus_id']);c=gate[r['locus_id']];desc='nr 中病毒报告恢复' if a['contribution']=='viral_reporting_recovery_in_uniform_DIAMOND' else '病毒专库带来的补充检出' if a['contribution']=='viral_database_panel_increment_in_uniform_DIAMOND' else '非病毒报告独立贡献待解释';lines.append(f"| {r['accession']} | {r['species']} | {desc} | {c.get('tree_primary_class','')} | 合格远缘同源中病毒优先；方向未证明 |")
lines+=['','“合格远缘同源”排除 Self/Close，要求 identity≥30%、完整查询覆盖≥50%、E-value≤1e-10。病毒最高同源、病毒作为载体、历史转移方向是不同问题。未发现合格细胞远缘参照也不等于没有任何细胞生物同源。','','## 从 19 条记录到最终集合','','- 19 条记录对应 17 个选定的注释蛋白位置、16 个不同蛋白 accession。两个核酸 ORF 与注释蛋白同区域、同方向、同读框，在新 DIAMOND 结果出来前已固定为支持记录，不作为独立新增，也不在蛋白版本失败后换版本补过关。','- 统一 DIAMOND 初筛通过 6 条记录：5 条注释蛋白和 1 条上述同位点核酸记录。EDR27775.1 不通过原首个定位 contig 的要求，余下 4 条注释蛋白进入最终树和去重复核。',f'- 最终新增 {n} 条。逐条判据与未纳入原因见下表和 [完整流程表](tables/final_candidate_gate19.tsv)。','','| 原始记录／蛋白 | 最终状态 | 未通过或未执行的阶段 |','|---|---|---|']
for z in gates:
 failed=[k+': '+v for k,v in z.items() if isinstance(v,str) and v.startswith('FAIL')]
 reason='；'.join(failed)
 if z.get('initial_failure_reasons'):reason+='（初筛：'+z['initial_failure_reasons']+'）'
 if z['final_status']=='supporting_alternate_ORF_not_independent_addition':reason='同位点支持版本，不独立计数；原核酸链未独立认证'
 lines.append(f"| {z['locus_id']} / {z['accession']} | {z['final_status']} | {reason or ('全部所需阶段通过' if z['all_required_stages_passed']=='True' else z['nonpass_stages'])} |")
lines+=['','这里“未纳入”包括固定门槛不通过、重复支持记录或未闭合证据，不等于证明该基因不曾发生 HGT。早期明确失败的记录没有继续当作合格对象建树。','','## 方法与边界','','### 统一 DIAMOND 搜索与初筛','','按本次用户确认，采用 DIAMOND，而不是声称完全重现原 MMseqs 搜索。19 个固定查询分别搜索同一本地 nr 去冗余蛋白库与同一合并病毒库：','','```text','diamond blastp --sensitive --evalue 1e-10 --max-target-seqs 2000','```','','搜索阶段不设置 identity 或 query-cover 下限，保留弱近缘参照；nr 用 24 线程，病毒库用 2 线程。完整命令、数据库大小和查询 SHA-256 见 [搜索记录](diamond_review/search_complete.json)；另有 [软件版本](diamond_review/software_version.txt) 和 [分析规则](diamond_review/analysis_plan.json)。这些查询是对现有候选的复核，不是重新对全部 105 物种开展一轮 de novo 搜索。','','初筛沿用原脚本的 Out 比例、数量、Self、AI、第一 Out 的 identity 和 coverage 条件。AI 为第一 Out 与第一 Close 的 identity×alignment length 之差；原 coverage 定义为第一 Out 比对长度除以第一 Self 比对长度，不能与完整 query coverage 混称。普通分支要求 Out_pct>0.9、Out>1、Out≥Self、Self>0、AI>0、第一 Out identity≥0.3、coverage≥0.7，并通过熵门槛；细菌特殊分支按原代码 coverage≥0.5，可绕过熵条件。','','**熵阈值使用此前全物种背景计算并冻结的、病毒分类修复后各物种阈值，作为外部校准。没有用 19 条候选重估阈值，也没有重新对全部蛋白组执行 DIAMOND 后校准。** 因而本次是固定历史校准下的补充筛选，不能表述为对整套 DIAMOND 背景重新估计了完全相同的 KDE 分布。','','对同一批 DIAMOND 查询比较三条分支：nr 中不计病毒、nr 中计入病毒、nr＋病毒专库。它们的初筛通过数分别为 2、5、6（含核酸重复记录）；最终纳入的增量归属在 [对照分支表](tables/accepted_candidate_attribution.tsv) 中逐条给出。增加数据库不会保证每个旧阳性都仍通过，因为第一 Out、分类分布和熵都可能改变。','','XP_004257338.1 和 CAG4718425.1 达到 2,000 hit 上限，均另有第一 Out coverage 等明确门槛失败。本轮不认证这两条，也不把它们截断后的计数解释为完整同源物种分布。若未来重新争取纳入，需要扩大命中上限并重新核查。','','### 后续过滤、取样与系统树','','沿用已复核的五个 BUSCO 谱系检查、EGT/共生来源标记、原基因组覆盖及首个定位 contig 规则、修订阶段 contig tier 和双侧真核侧翼证据。动物保守过滤以新第一 Out 的 identity/coverage 重算比较，不能沿用旧第一 Out 的差值。所有 19 条已保存基因组区域的坐标和方向均曾精确复现；这不等价于证明唯一拷贝或每条 GFF 模型完全正确。','','树取样复用原动态配额及属／种多样性选择函数，保留长度≥80 aa、identity≥20% 的参照；查询和同源片段按原规则加入 20 aa 侧翼、去除相同序列。17 条蛋白查询全部与旧搜索输入全序列一致，因此旧门槛下合格参照作为独立来源标签补入树背景，不混入新的 DIAMOND 初筛计数。这避免仅因新搜索过滤而丢失已知 Close 参照。','','MAFFT、trimAl automated1 后，原模型 LG+F+G4 与固定 LG+C20+F+G4 敏感性模型均使用 1,000 UFBoot 和 1,000 SH-aLRT，随机种子 20260915。原分类器及类别白名单保持不变，并检查收敛、长枝、查询/树叶/比对一致性。IQ-TREE 将标签内等号改为下划线的处理已验证为一一对应，不属于生物学阈值放宽。最后一棵 C20 从检查点增加线程续跑，资源变化另有记录。','','“未触发原长枝判据”不能写成已绝对排除长枝吸引。保留无 Close、局部低支持、混合邻域及混合模型过拟合警告；两种模型并不构成彼此独立的实验验证。新增仍按原允许的候选级别报告，不升级成全部强证据 HGT。','','另核查病毒与细胞 Out 的树证据，排除 Self/Close 对供体距离的干扰。若两模型均给出强支持的细胞 Out 侧分布并明显更近，则不按病毒优先来源纳入；其他情形仍仅称潜在病毒来源。具体值见 [来源树检查](tables/viral_source_topology_review.tsv)。','','### 去重与旧目录保留','','90% CD-HIT 使用声明的“冻结历史分类修复初筛池＋本轮初筛通过蛋白”输入，并保留标题精确对应及缺失检查；这不是重新计算全部 105 蛋白组的 DIAMOND 初筛池。后续使用原物种内 `mmseqs easy-cluster --min-seq-id 0.8 -c 0.8 --cov-mode 1`，对历史预聚类池与最终合格新候选共同处理。再检查与原 8,985 的代表关系及基因组坐标重叠。两个不同 contig 上的同一 EAX69617.1 查询仍保留双位置记录，但没有据此宣称两个真实独立拷贝。','','## Irwin 比较','','按作者展开的最终事件成员及 Table 1 Direction 计数，主比较保持污染分数≥2。全长一致指完整氨基酸序列相等，对双方都是全覆盖、无氨基酸差异；不只是长度相等，也不要求 accession 相同。','','| 集合 | 作者方向 | 不同 Irwin 蛋白 | 不同作者事件 |','|---|---|---:|---:|']
for r in overlap:lines.append(f"| {r['cohort']} | {r['direction']} | {r['unique_Irwin_proteins']} | {r['unique_Irwin_events']} |")
lines+=['','S2 中的序列相同但没有明确最终事件成员关系，不能自动归为 EtoV，也不能作为作者最终 HGT 的重叠证据。方向和事件的逐条支持见 [Irwin 证据表](tables/Irwin_final_overlap_proof.tsv)。','','## 交付文件','',f'- [合并目录：{len(combined):,} 条](tables/{name})：原 8,985 条及其原字段完全保留；新增带本轮方法和来源解释标记。','- [新增完整特征与审查信息](tables/accepted_new_candidates.tsv)：保留此前 GC、CDS、最低年龄、Ka/Ks 等字段，并增加本轮结果。此前 Ka/Ks 的供体若与本轮最佳同源不同，单列一致性标记，不冒充已针对新供体重算。','- [新增蛋白 FASTA](accepted_new.proteins.faa) 与 [对应方向的基因组区域](accepted_new.regions.fna)。','- [最终逐条判定](tables/final_candidate_gate19.tsv)、[结构化结果](final_manifest.json)、[统一初筛三分支](diamond_review/initial_screen_three_arms.tsv)。','','本轮没有生物安全拒绝或因此跳过的任务。两条核酸支持记录没有被宣称独立通过完整原核酸筛选链；它们也没有作为新增条目计数。']

insert=lines.index('## 交付文件')
detail=['## 证据强度及使用注意','',
'4 条通过的是本轮声明的候选纳入条件，并非 4 条均已确定病毒→寄生虫的方向。三条 Acanthamoeba 候选在本轮合格远缘命中中没有细胞参照；较宽的树参照中，XP_004339817.1 仍有 1 条细胞 Out。RNC50057.1 则是病毒最佳合格远缘同源，但原模型中最近细胞 Out 距离略短于病毒，C20 中次序相反，保留来源不确定。','',
'四条均无 Close 组；原树类别分别为 3 条 possible_HGT_no_close_donor_leaning 和 1 条 recipient_lineage_like_no_close_possible_ancestral_HGT。ELR15564.1 保留 query_branch_high_ratio，XP_004339817.1 保留 local_split_low_support；四条 C20 都出现部分混合权重接近零的过拟合警告。它们未触发原长枝排除规则，不能写成已彻底排除长枝吸引。','',
'新增中 2 条属于统一 DIAMOND 下的 nr 病毒报告恢复，2 条需要病毒专库才通过本轮初筛。这是本次三分支比较的归因，不能据此声称原 MMseqs 历史搜索的唯一差异就是增加数据库。','',
'旧目录的病毒最高同源数量沿用已完成的 74 条统一复核：50 条病毒更高、1 条无合格细胞远缘参照，其中 7 条另有可竞争的 UNKNOWN。加上本轮 4 条，得到这一已复核子集内 55 条病毒优先／无合格细胞远缘参照记录；不能将 55 写成全 8,989 条均被同参数重搜后得到的完整总数。原 1,194 条病毒阳性注释仍保留。','',
'合并表保留旧 Irwin 字段供追溯；当前明确事件方向应读取 final_review_Irwin_* 新字段，不能再以旧的 missing_expanded_table 状态作结论。','']
lines[insert:insert]=detail
(F/'REPORT_zh.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(new=n,total=len(combined),final=g['final_count_determined_within_declared_scope'],catalogue=name,Irwin=overlap),indent=2))

# Apply the subsequent support-stratified Irwin correction on every regeneration.
import subprocess,sys
correction=F/"irwin_support_audit_20260915/scripts/publish_correction.py"
if correction.exists():
 subprocess.run([sys.executable,str(correction)],check=True)
