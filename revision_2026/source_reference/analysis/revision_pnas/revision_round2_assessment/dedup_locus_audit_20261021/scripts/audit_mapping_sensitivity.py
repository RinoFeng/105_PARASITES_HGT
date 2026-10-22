from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,json,re
import numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent;G=R/'innovation_execution/genome_integration';T=P/'tables'
csv.field_size_limit(10000000)
def rows(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:yield from csv.DictReader(f,delimiter='\t')
def toks(x):return {z for z in re.split('[;,|]',x or '') if z and z not in {'NA','nan','-'}}
def yes(x):return str(x).lower() in {'true','1'}
def write(n,rr):
 rr=list(rr)
 with (T/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]) if rr else ['none'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
panel={r['species']:r for r in rows(T/'09_species_counts_and_sensitivity.tsv')};inv={r['unit_id']:r for r in rows(R/'manuscript_evidence_completion_20261003/data/final48_locus_inventory.tsv.gz')};final={r['hgt_id'] for r in rows(R/'innovation_final8989_20260915/inputs/final_catalogue_8989.tsv')}
amb=[r for r in rows(G/'refinement_20260914/candidate_locus_evidence_audit.tsv') if r['record_id'] in final and r['species'] in panel and r['direct_evidence_status']=='multiple_own_record_gene_links'];need={u for r in amb for u in toks(r['own_record_gene_ids'])};og={r['unit_id']:r['inferred_family'] for r in rows(G/'refinement_20260914/unit_orthogroup_membership.tsv.gz')};models={}
for inp in rows(G/'species_inputs.tsv'):
 if inp['species'] not in {r['species'] for r in amb}:continue
 p=G/'measurements'/(inp['species_id']+'.tsv.gz')
 for r in rows(p):
  if r['unit_id'] in need:models[r['unit_id']]={k:r.get(k,'') for k in ['unit_id','species','sequence_id','gene_id','contig','strand','coding_start_0','coding_end_0','coding_structure','eligible_structure','direct_catalogue_records']};models[r['unit_id']]['OG']=og.get(r['unit_id'],'')
ambrows=[];add={}
for r in amb:
 uu=toks(r['own_record_gene_ids']);mm=[models[u] for u in uu if u in models];ogset={x['OG'] for x in mm};seqs={x['sequence_id'] for x in mm};pairnonoverlap=all(a['contig']!=b['contig'] or int(a['coding_end_0'])<=int(b['coding_start_0']) or int(b['coding_end_0'])<=int(a['coding_start_0']) for i,a in enumerate(mm) for b in mm[i+1:]);ok=len(mm)==len(uu) and all(yes(x['eligible_structure']) for x in mm) and len(ogset)==1 and '' not in ogset and pairnonoverlap
 if ok:
  for x in mm:
   if x['unit_id'] not in inv:add[x['unit_id']]=dict(x,anchor_record=r['record_id'],HGT_confirmed_by_this_audit=False)
 ambrows.append(dict(record_id=r['record_id'],species=r['species'],clade=panel[r['species']]['clade'],gene_models=len(uu),models_recovered=len(mm),eligible_models=sum(yes(x['eligible_structure']) for x in mm),distinct_sequences=len(seqs),distinct_OGs=len(ogset),OGs=';'.join(sorted(ogset)),pairwise_nonoverlapping=pairnonoverlap,conditional_addition_eligible=ok,units=';'.join(sorted(uu))))
write('12_multiple_locus_candidate_records.tsv',ambrows);write('13_multiple_locus_conditional_units_NOT_NEW_HGT.tsv',add.values())
# Apply the identical-sequence unique-locus criterion symmetrically to the background.
byseq=defaultdict(list)
for u,r in inv.items():byseq[r['species'],r['sequence_id']].append(u)
rm={u for uu in byseq.values() if len(uu)>1 for u in uu if inv[u]['role']=='background'}
# Counts and comparisons retain exactly the same 48 species and frozen 1,000 draws.
extra90={r['unit_id']:r for r in rows(T/'07_potential_uncounted_loci_NOT_NEW_HGT.tsv')}
scenarios={'published':(set(),{}),'background_unique_sequence_locus_only':(rm,{}),'conditional_identical_multi_locus_candidates':(set(),add),'conditional90_and_identical_multi_locus':(set(),dict(extra90,**add))}
spout=[];contr=[];clades=['Amoebozoa','Apicomplexa','Euglenozoa','Metazoa']
for name,(removed,extras) in scenarios.items():
 vals=[]
 for sp in sorted(panel):
  c=[dict(x) for u,x in inv.items() if x['species']==sp and x['role']=='candidate'];b=[x for u,x in inv.items() if x['species']==sp and x['role']=='background' and u not in removed];c.extend(x for u,x in extras.items() if x['species']==sp)
  rr=dict(scenario=name,species=sp,clade=panel[sp]['clade'],C=len({x['OG'] for x in c}),N=len(c),B=len({x['OG'] for x in b}),M=len(b));vals.append(rr);spout.append(rr)
 rng=np.random.default_rng(20260922);arr={}
 for cl in clades:
  z=[r for r in vals if r['clade']==cl];aa=[]
  for k in range(1000):
   zz=[z[i] for i in rng.choice(len(z),3,replace=False)];c,b,n,m=[sum(r[j] for r in zz) for j in ['C','B','N','M']];aa.append((c/b,(n/c)/(m/b),n/m))
  arr[cl]=np.asarray(aa)
 for cl in clades[:-1]:
  rat=arr[cl]/arr['Metazoa']
  for j,key in enumerate(['family_ratio','within_family_loci_ratio','locus_ratio']):
   v=rat[:,j];contr.append(dict(scenario=name,clade=cl,metric=key,fold=float(np.exp(np.log(v).mean())),q025=float(np.quantile(v,.025)),q975=float(np.quantile(v,.975)),formal_P_value=False))
write('14_mapping_sensitivity_species.tsv',spout);write('15_mapping_sensitivity_contrasts.tsv',contr)
# 50% cluster versus gene family partition: evaluate only strict linked loci.
fifty=list(rows(T/'11_strict_panel_50cluster_OG_links.tsv'));cl2og=defaultdict(set);og2cl=defaultdict(set);cl2u=defaultdict(set)
for r in fifty:cl2og[r['cluster50_rep']].add(r['OG']);og2cl[r['OG']].add(r['cluster50_rep']);cl2u[r['cluster50_rep']].add(r['unit_id'])
write('16_50cluster_vs_orthogroup.tsv',[dict(cluster50_rep=c,orthogroups=len(v),OGs=';'.join(sorted(v)),candidate_loci=len(cl2u[c])) for c,v in sorted(cl2og.items())])
summary=dict(multiple_mapping_records=len(amb),multiple_mapping_models=len(need),multiple_mapping_eligible_records=sum(r['conditional_addition_eligible'] for r in ambrows),conditional_additional_models=len(add),ambiguous_by_clade=dict(Counter(r['clade'] for r in ambrows)),conditional_additional_models_by_clade=dict(Counter(panel[r['species']]['clade'] for r in add.values())),background_removed_singleton_sensitivity=len(rm),background_removed_by_clade=dict(Counter(panel[inv[u]['species']]['clade'] for u in rm)),global50_strict_loci=len({r['unit_id'] for r in fifty}),global50_strict_nonrep_loci=len({r['unit_id'] for r in fifty if r['is_rep']=='False'}),global50_clusters=len(cl2og),global50_clusters_spanning_multiple_OGs=sum(len(v)>1 for v in cl2og.values()),OGs_in50_map=len(og2cl),OGs_spanning_multiple50_clusters=sum(len(v)>1 for v in og2cl.values()),strict_linked_records_data_types=dict(Counter('nucleotide_coordinate' if ':' in r['accession'] else 'protein_accession' for r in fifty)))
(P/'mapping_sensitivity_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
