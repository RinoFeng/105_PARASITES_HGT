"""Fixed-scope DIAMOND recheck using frozen whole-species entropy calibration.
Preserve all raw hits. Same raw query and parameters in nr/panel arms.
"""
from pathlib import Path
from collections import defaultdict,Counter
import csv,json,sys,math,subprocess,hashlib
F=Path(__file__).resolve().parents[1];D=F.parent;B=D.parent;S=F/'diamond_review'
assert json.load((S/'search_complete.json').open())['complete']
sys.path.insert(0,str(B/'scripts'));import merge_catalog_viral_ranking as M
from viral_taxonomy import ViralTaxonomy
M.D=S
records=M.rows(D/'tables/viral_priority19_all_evidence.tsv');record={r['locus_id']:r for r in records};spec={r['genome']:r for r in M.rows(B/'inputs/species_manifest.tsv')}
raw=[];nrids=set();vids=set();counts=Counter()
for panel in ['nr','viral']:
 for z in csv.reader((S/(panel+'.tsv')).open(),delimiter='\t'):
  raw.append((panel,z));(nrids if panel=='nr' else vids).add(z[1]);counts[panel,z[0]]+=1
ids=S/'nr_ids.txt';ids.write_text('\n'.join(sorted(nrids))+'\n');taxfile=S/'nr_taxonomy.tsv'
if not taxfile.exists():
 subprocess.run([str(B/'reports/merged_viral_ranking/filter_saved_query_rows'),str(ids),'EXTERNAL_DATA_ROOT/Data/DB/nr/nr_cluster_taxid_formatted_final.tsv',str(taxfile)],check=True)
nrmap={z[0]:dict(taxids=z[1],lca=z[2],detail=z[3]) for z in csv.reader(taxfile.open(),delimiter='\t') if len(z)>=4}
needed={spec[r['genome']]['taxid'] for r in records}
for x in nrmap.values():needed.update(x['taxids'].split(';'));needed.add(x['lca'])
tax=M.taxonomy_load(needed);vt=ViralTaxonomy(B,vids,tax);vt.enrich(B,set().union(*vt.taxids.values()) if vt.taxids else set())
(S/'resolved_taxonomy.json').write_text(json.dumps(tax))
hits=defaultdict(list)
for panel,z in raw:
 q,acc=z[:2];tids=vt.taxids.get(acc,set()) if panel=='viral' else set(nrmap.get(acc,{}).get('taxids','').split(';'))-{''};tid=vt.lca.get(acc,'') if panel=='viral' else nrmap.get(acc,{}).get('lca','');sid=spec[record[q]['genome']]['taxid'];groups={M.get_group(t,sid,tax) for t in tids};group='VIRAL' if panel=='viral' else next(iter(groups)) if len(groups)==1 else 'UNKNOWN'
 h=dict(locus_id=q,subject=acc,source=panel,taxid=tid,all_source_taxids=';'.join(sorted(tids)),group=group,lineage=tax.get(tid,{}).get('detail',''),taxonomy_short=tax.get(tid,{}).get('short',''),bitscore=float(z[9]),pident=float(z[2]),alignment_length=int(z[3]),qstart=int(z[4]),qend=int(z[5]),sstart=int(z[6]),send=int(z[7]),evalue=float(z[8]),query_length=int(z[10]),subject_length=int(z[11]),query_coverage=(abs(int(z[5])-int(z[4]))+1)/int(z[10]))
 hits[q].append(h)
M.write(S/'all_hits_annotated.tsv',[h for hh in hits.values() for h in hh])
results=[]
for r in records:
 q=r['locus_id'];cal=json.load((B/'recovery/rounded_legacy/protein'/(r['genome']+'.summary.json')).open());cut=cal['fixed_entropy_cutoff']
 for arm in ['nr_viral_unreported','nr_viral_reported','nr_plus_viral_panel']:
  hh=[h for h in hits[q] if arm=='nr_plus_viral_panel' or h['source']=='nr'];hh=sorted(hh,key=lambda h:(-h['bitscore'],-h['alignment_length'],-h['pident'],h['subject'],h['source']));unique={}
  for h in hh:
   key=(h['subject'],h['taxid']);unique.setdefault(key,h)
  hh=list(unique.values());taxhits={}
  for h in hh:taxhits.setdefault(h['taxid'] or 'UNRESOLVED_'+h['group'],h)
  known=list(taxhits.values());n=Counter();classes=Counter();first={}
  for h in known:
   group=h['group'];short=h['taxonomy_short'].split(';')
   if len(short)>2:classes[short[2]]+=1
   if group=='VIRAL' and arm=='nr_viral_unreported':continue
   if group=='UNKNOWN':n['UNKNOWN']+=1;continue
   if group in ['VIRAL','CELLULAR_OUT']:group='OUT'
   n[group]+=1;first.setdefault(group,h)
  out=n['OUT'];close=n['CLOSE'];selfn=n['SELF'];ratio=round(out/(out+close),3) if out+close else 0.;ct=sum(classes.values());entropy=-sum(v/ct*math.log2(v/ct) for v in classes.values()) if ct else 0.;o=first.get('OUT');c=first.get('CLOSE');s=first.get('SELF');coverage=o['alignment_length']/s['alignment_length'] if o and s else 0.;identity=o['pident']/100 if o else 0.;ai=(o['pident']/100*o['alignment_length'] if o else 0)-(c['pident']/100*c['alignment_length'] if c else 0)
  dist=ratio>.9 and out>1 and out>=selfn and selfn>0 and ai>0 and identity>=.3
  normal=dist and coverage>=.7 and cut is not None and entropy<cut
  special=dist and coverage>=.5 and bool(o and 'Bacteria' in o['lineage'])
  qual=[h for h in hh if h['pident']>=30 and h['query_coverage']>=.5]
  v=next((h for h in qual if h['group']=='VIRAL'),None);cell=next((h for h in qual if h['group']=='CELLULAR_OUT'),None);unk=next((h for h in qual if h['group']=='UNKNOWN'),None)
  rank=M.status(v['bitscore'] if v else None,cell['bitscore'] if cell else None)
  failures=[]
  for ok,label in [(ratio>.9,'Out_pct<=0.9'),(out>1,'Out_n<=1'),(out>=selfn,'Out_n<Self_n'),(selfn>0,'no_SELF'),(ai>0,'AI<=0'),(identity>=.3,'first_Out_identity<0.3'),(coverage>=.5 if special else coverage>=.7,'first_Out_coverage'),(special or cut is not None and entropy<cut,'entropy')]:
   if not ok:failures.append(label)
  z=dict(locus_id=q,accession=r['accession'],genome=r['genome'],workflow=r['workflow'],arm=arm,raw_deduplicated_hits=len(hh),unique_taxon_buckets=len(known),Self_num=selfn,Close_num=close,Out_num=out,Unknown_num=n['UNKNOWN'],Out_pct=ratio,first_Out_subject=o['subject'] if o else '',first_Out_group=o['group'] if o else '',first_Out_coverage=coverage,first_Out_identity=identity,AI=ai,entropy=entropy,frozen_species_entropy_cutoff=cut,initial_distribution_pass=normal or special,initial_failure_reasons=';'.join(failures),viral_rank_qualified=rank,best_qualified_viral_subject=v['subject'] if v else '',best_qualified_viral_bitscore=v['bitscore'] if v else '',best_qualified_cellular_subject=cell['subject'] if cell else '',best_qualified_cellular_bitscore=cell['bitscore'] if cell else '',unknown_competes_with_best_distant=bool(unk and (v or cell) and unk['bitscore']>=max(v['bitscore'] if v else 0,cell['bitscore'] if cell else 0)),nr_hit_cap_reached=counts['nr',q]>=2000,viral_hit_cap_reached=counts['viral',q]>=2000,calibration='frozen_external_legacy_whole_species_viral_repaired_threshold;not_reestimated_from19',nucleotide_note='translated_ORF_recheck_not_complete_original_nt_window_chain' if r['workflow']=='nucleotide' else '')
  results.append(z)
M.write(S/'initial_screen_three_arms.tsv',results)
combined=[r for r in results if r['arm']=='nr_plus_viral_panel'];M.write(S/'initial_screen_combined19.tsv',combined)
summary=dict(records=19,arm_pass_counts={arm:sum(r['initial_distribution_pass'] for r in results if r['arm']==arm) for arm in ['nr_viral_unreported','nr_viral_reported','nr_plus_viral_panel']},hit_cap_records=[r['locus_id'] for r in combined if r['nr_hit_cap_reached'] or r['viral_hit_cap_reached']],calibration_limitation='fixed external historical species threshold; no new whole-proteome recalibration',full_nt_chain_completed=False)
(S/'analysis_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
