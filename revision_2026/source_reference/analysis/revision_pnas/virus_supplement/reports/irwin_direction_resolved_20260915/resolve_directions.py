from pathlib import Path
from collections import defaultdict,Counter
import pandas as pd,numpy as np,json,hashlib,math
R=Path(__file__).resolve().parent;B=R.parent.parent;S=R/'source_tables'
final=pd.read_csv(S/'TableS1_HGTSummary_expanded.tab',sep='\t',keep_default_na=False,dtype=str)
raw=pd.read_csv(S/'HGT_summary_table_noContaminationCleaning.tab',sep='\t',keep_default_na=False,dtype=str)
x=pd.read_excel(B/'literature/irwin2022/supplementary_tables_1_5.xlsx',sheet_name=0,header=6,keep_default_na=False)
p=pd.read_csv(B/'reports/irwin_direction_audit/S2_all1423_direction_join_pending.tsv',sep='\t',keep_default_na=False)
trace=pd.read_csv(B/'reports/irwin_913_filter_trace/literature_protein_stage_trace_913.tsv',sep='\t',keep_default_na=False)
pidx=p.set_index('irwin_table2_id').to_dict('index');target=set(pidx)
parents={};ranks={};merged={}
with (B/'databases/taxonomy/nodes.dmp').open() as f:
 for l in f:
  a=l.split('|');parents[a[0].strip()]=a[1].strip();ranks[a[0].strip()]=a[2].strip()
with (B/'databases/taxonomy/merged.dmp').open() as f:
 for l in f:
  a=l.split('|');merged[a[0].strip()]=a[1].strip()
cache={}
def sp(t):
 if t in cache:return cache[t]
 orig=t;seen=set()
 while t in merged and t not in seen:seen.add(t);t=merged[t]
 seen=set()
 while t and t not in seen:
  if ranks.get(t)=='species':cache[orig]=t;return t
  seen.add(t);t=parents.get(t,'')
 cache[orig]='';return ''
m=pd.read_csv(B/'inputs/species_manifest.tsv',sep='\t',dtype=str,keep_default_na=False);study={sp(r.taxid) for r in m.itertuples()};study.discard('')
def ids(v):return {x.strip() for x in str(v).split(',') if x.strip().split('.',1)[0].isdigit() and '.' in x}
def isstudy(s):return sp(s.split('.',1)[0]) in study
def key(r,cols):return tuple(str(r[c]).strip().lower() for c in cols)
xg=defaultdict(list)
for i,r in x.iterrows():xg[key(r,['Cluster','Direction','Recipient_taxaID','Donor_taxaID'])].append((i+8,r))
# Stable events use full source row content; exact duplicate rows are retained as source aliases.
groups={}
for i,r in final.iterrows():
 h=hashlib.sha256('\t'.join(r).encode()).hexdigest()[:20]
 if h not in groups:groups[h]={'record':r,'rows':[]}
 groups[h]['rows'].append(i+2)
rawg=defaultdict(list)
for i,r in raw.iterrows():rawg[(r.Cluster,r.direction)].append((i+2,r,ids(r.recipients)))
ev=[];links=[];beforelinks=[];beforematches=[]
for h,g in groups.items():
 r=g['record'];eid='Irwin_'+h;rec=ids(r.recipients);don=ids(r.donors);side=don if r.direction=='EtoV' else rec
 cand=xg[key(r,['Cluster','direction','taxid_recip','taxid_donor'])]
 match=[(row,a) for row,a in cand if math.isclose(float(a.Recipient_branch_length),float(r.branch_length),rel_tol=1e-10,abs_tol=1e-9)]
 assert match,(eid,r.Cluster)
 s1rows=sorted(row for row,a in match)
 # Any unresolved S1 ambiguity remains an alias list; the directions and event members come from expanded table itself.
 rawmatch=[(row,a) for row,a,rr in rawg[(r.Cluster,r.direction)] if rec.issubset(rr) and math.isclose(float(a.branch_length),float(r.branch_length),rel_tol=1e-10,abs_tol=1e-9)]
 for row,a in rawmatch:
  beforematches.append(dict(event_id=eid,final_expanded_rows=';'.join(map(str,g['rows'])),raw_row=row,raw_recipient_members=len(ids(a.recipients)),final_recipient_members=len(rec),removed_recipient_members=';'.join(sorted(ids(a.recipients)-rec)),mapping_rule='same_cluster_direction_branch_and_final_recipient_subset'))
 conf='weakly_supported' if r.taxid_donor=='Weakly supported' or r.name_donor=='Weakly supported' else 'not_marked_weak_in_expanded'
 studyids=sorted(s for s in side if isstudy(s));s2ids=sorted(side&target)
 row=dict(event_id=eid,expanded_rows=';'.join(map(str,g['rows'])),source_row_count=len(g['rows']),S1_excel_rows=';'.join(map(str,s1rows)),S1_mapping_unique_except_exact_duplicate=len(s1rows)==len(g['rows']),Cluster=r.Cluster,direction=r.direction,taxid_recip=r.taxid_recip,name_recip=r.name_recip,taxid_donor=r.taxid_donor,name_donor=r.name_donor,notes=r.notes,confidence=conf,support1=r.support1,support2=r.support2,recipient_sequence_count=len(rec),donor_sequence_count=len(don),eukaryotic_side_members_available=bool(side),study_species_member_count=len(studyids),S2_study_member_count=len(s2ids),raw_matches=len(rawmatch),member_role='eukaryotic_donor' if r.direction=='EtoV' else 'eukaryotic_side_direction_unknown' if r.direction=='unkn' else 'eukaryotic_recipient')
 ev.append(row)
 for sid in sorted(set(studyids)|set(s2ids)):
  links.append(dict(event_id=eid,Cluster=r.Cluster,direction=r.direction,confidence=conf,member_role=row['member_role'],author_sequence_id=sid,in_original_S2_1423=sid in target,taxonomy_in_105=isstudy(sid)))
for i,r in raw.iterrows():
 for sid in ids(r.donors if r.direction=='EtoV' else r.recipients)&target:beforelinks.append(dict(raw_row=i+2,Cluster=r.Cluster,direction=r.direction,irwin_table2_id=sid))
e=pd.DataFrame(ev);l=pd.DataFrame(links);br=pd.DataFrame(beforelinks);bm=pd.DataFrame(beforematches)
e.to_csv(R/'final_events_with_S1_correspondence.tsv',sep='\t',index=False);bm.to_csv(R/'final_to_prefilter_event_correspondence.tsv',sep='\t',index=False);br.to_csv(R/'prefilter_event_S2_members.tsv',sep='\t',index=False)
assert len(e)==6740 and e.source_row_count.sum()==6741
assert e.S1_mapping_unique_except_exact_duplicate.all()
joined=l.merge(p,left_on='author_sequence_id',right_on='irwin_table2_id',how='left',validate='many_to_one',suffixes=('','_old_pending'))
joined=joined.drop(columns=['direction_old_pending'],errors='ignore');joined.to_csv(R/'final_event_study_protein_relations.tsv',sep='\t',index=False)
s2=l[l.in_original_S2_1423].merge(p,left_on='author_sequence_id',right_on='irwin_table2_id',validate='many_to_one',suffixes=('','_old_pending')).drop(columns=['direction_old_pending'])
s2.to_csv(R/'final_event_S2_protein_relations.tsv',sep='\t',index=False)
mp=defaultdict(list)
for rr in s2.to_dict('records'):mp[rr['irwin_table2_id']].append(rr)
rawmap=br.groupby('irwin_table2_id').direction.agg(lambda s:';'.join(sorted(set(s)))).to_dict()
pr=[]
for rr in p.to_dict('records'):
 rel=mp[rr['irwin_table2_id']];dirs=sorted({a['direction'] for a in rel});rr['direction']=';'.join(dirs);rr['Irwin_final_event_membership']='explicit_final_event_member' if rel else 'not_listed_on_eukaryotic_side_in_final_expanded_table';rr['final_event_ids']=';'.join(sorted({a['event_id'] for a in rel}));rr['raw_directions']=rawmap.get(rr['irwin_table2_id'],'');rr['in_final_VtoE']='VtoE' in dirs;rr['in_final_unkn']='unkn' in dirs;rr['in_final_EtoV']='EtoV' in dirs;rr['final_membership_status']='listed_final' if rel else 'listed_prefilter_only' if rr['raw_directions'] else 'absent_from_both_event_member_lists';pr.append(rr)
pr=pd.DataFrame(pr);pr.to_csv(R/'S2_1423_direction_resolved.tsv',sep='\t',index=False)
hit=pr[pr.original_catalogue_exact|pr.supplement_new_exact].merge(trace[['irwin_table2_id','species','genome','accession','new_retained_locus_ids']],on='irwin_table2_id',validate='one_to_one');hit.to_csv(R/'42_original_overlap_direction_audit.tsv',sep='\t',index=False)
summary=[];filters=[];eventrows=[]
for direction in ['VtoE','unkn','EtoV']:
 for confidence in ['all','not_marked_weak_in_expanded','weakly_supported']:
  a=l[(l.direction==direction)&((l.confidence==confidence) if confidence!='all' else True)];b=s2[(s2.direction==direction)&((s2.confidence==confidence) if confidence!='all' else True)];q=pr[pr.irwin_table2_id.isin(b.irwin_table2_id)];exact=q[q.exact_in_local_proteome];old=q[q.original_catalogue_exact];new=q[q.supplement_new_exact];union=q[q.original_catalogue_exact|q.supplement_new_exact]
  total_events=set(a.event_id);s2events=set(b.event_id);covered=set(b.loc[b.original_catalogue_exact|b.supplement_new_exact,'event_id']);comparable=set(b.loc[b.exact_in_local_proteome,'event_id']);old_events=set(b.loc[b.original_catalogue_exact,'event_id']);new_events=set(b.loc[b.supplement_new_exact,'event_id'])
  summary.append(dict(direction=direction,confidence=confidence,events_with_105_members=len(total_events),unique_105_member_sequences=a.author_sequence_id.nunique(),events_with_S2_1423_members=len(s2events),S2_unique_proteins=len(q),S2_exact_local_proteins=len(exact),S2_old_overlap_proteins=len(old),S2_new_overlap_proteins=len(new),S2_union_overlap_proteins=len(union),old_covered_events=len(old_events),new_covered_events=len(new_events),union_covered_events=len(covered),events_with_exact_comparable_members=len(comparable),exact_protein_match_fraction=len(union)/len(exact) if len(exact) else None,overlap_with_tree_flags=int(union.latest_tree_outside_original_whitelist.sum())))
  for stagecol in ['availability_state','historical_terminal_stage','new_terminal_stage']:
   use=q if stagecol=='availability_state' else exact
   for stage,n in use[stagecol].value_counts().items():filters.append(dict(direction=direction,confidence=confidence,stage_type=stagecol,stage=stage,unique_proteins=int(n),scope='S2_1423_final_event_members'))
for eid,a in s2.groupby('event_id'):
 q=a.drop_duplicates('irwin_table2_id');ex=q[q.exact_in_local_proteome];covered=q[q.original_catalogue_exact|q.supplement_new_exact]
 eventrows.append(dict(event_id=eid,direction=q.direction.iloc[0],confidence=q.confidence.iloc[0],related_S2_proteins=len(q),exact_comparable_proteins=len(ex),covered_proteins=len(covered),at_least_one_covered=bool(len(covered)),all_comparable_covered=bool(len(ex)) and len(covered)==len(ex),all_related_S2_covered=len(covered)==len(q)))
summ=pd.DataFrame(summary);summ.to_csv(R/'direction_overlap_summary.tsv',sep='\t',index=False);pd.DataFrame(filters).to_csv(R/'direction_stratified_filter_counts.tsv',sep='\t',index=False);pd.DataFrame(eventrows).to_csv(R/'per_event_coverage.tsv',sep='\t',index=False)
# Preserve exact first-failure detail plus non-exclusive failure fields, scoped to final author direction.
t=s2.merge(trace,on='irwin_table2_id',how='inner',suffixes=('','_trace'),validate='many_to_one');t.to_csv(R/'direction_event_protein_filter_trace.tsv.gz',sep='\t',index=False,compression='gzip')
stat={'complete_event_member_join':True,'expanded_source_rows':len(final),'unique_full_row_events':len(e),'exact_duplicate_source_rows':int(final.duplicated().sum()),'all_events_mapped_to_S1':True,'all_S1_mapping_unique_except_exact_duplicates':bool(e.S1_mapping_unique_except_exact_duplicate.all()),'final_events_with_prefilter_match':int((e.raw_matches>0).sum()),'final_events_without_prefilter_match':int((e.raw_matches==0).sum()),'42_overlap_direction_sets':hit.direction.value_counts().to_dict(),'42_overlap_membership_status':hit.final_membership_status.value_counts().to_dict(),'1423_membership_status':pr.final_membership_status.value_counts().to_dict(),'source_membership_not_replaced_by_family_or_tree_inference':True,'safety_skips':0}
(R/'join_status.json').write_text(json.dumps(stat,indent=2)+'\n');print(json.dumps(stat,indent=2));print(summ[summ.confidence=='all'].to_string(index=False))
