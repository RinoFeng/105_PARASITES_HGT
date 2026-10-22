"""Evaluate all five submitted Tier3 supports after both discovery routes are available."""
import csv,json,time,subprocess,shutil
from pathlib import Path
from collections import defaultdict
import filter_protein_discovery as P
from map_protein_candidates import tier
O=P.O;D=O/'reports'
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t')) if p.exists() and p.stat().st_size else []
def write(p,rr):
 with p.open('w') as f:
  if rr:
   keys=list(dict.fromkeys(k for r in rr for k in r));w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rr)
def yes(x):return str(x)=='True'
def overlaps(a,b):
 if a['strand']!=b['strand']:return False
 x,y=int(a['start_1based']),int(a['end_1based']);u,v=int(b['start_1based']),int(b['end_1based']);n=max(0,min(y,v)-max(x,u)+1)
 return 2*n>=y-x+1 and 2*n>=v-u+1
for rel in ['discovery/nucleotide_context/nucleotide_context.complete.json','discovery/missing_protein_genome_mapping/protein_mapping.complete.json','reports/iqtree_final_qc.complete.json','phylogeny/protein_recovered/protein_recovered_phylogeny.complete.json']:
 while not (O/rel).exists():time.sleep(15)
configs=[('protein','protein_context','protein_genome_mapping'),('protein_recovered','recovered_protein_context','recovered_protein_genome_mapping'),('protein_missing_context','missing_protein_context','missing_protein_genome_mapping')]
protein=[];byprotein=defaultdict(list);nucleotide=rows(O/'discovery/nucleotide_context/all_expanded_region_metrics.tsv');bynucleotide=defaultdict(list);best={}
for mode,context,mapping in configs:
 rr=rows(O/'discovery'/mapping/'all_protein_candidate_loci.tsv')
 for r in rr:
  r['_mode']=mode;r['_mapping']=mapping;protein.append(r)
  if yes(r['pass_mapping']):byprotein[(r['genome'],r['contig'])].append(r)
for r in nucleotide:
 if yes(r['pass_manuscript_coverage_variant']):bynucleotide[(r['genome'],r['contig'])].append(r)
for context in [c for m,c,d in configs]+['nucleotide_context']:
 for p in (O/'discovery'/context).glob('*.selected_homologs.tsv'):
  for h in rows(p):
   if h['group'] not in ['VIRAL','CELLULAR_OUT']:continue
   key=(h['genome'],h['query'])
   if key not in best or float(h['bitscore'])>float(best[key]['bitscore']):best[key]=h
changes=[]
for mode,rr,other in [('protein',protein,bynucleotide),('nucleotide',nucleotide,byprotein)]:
 for r in rr:
  old=yes(r['eligible_for_phylogeny']);q=r['accession'];h=best.get((r['genome'],q),{});lineage=h.get('taxonomy_detail') or P.T.get(h.get('taxid'),{}).get('detail','');ba='Bacteria;' in lineage or ';Archaea;' in lineage
  evidence=[x for x in other[(r['genome'],r['contig'])] if overlaps(r,x)];dual=bool(evidence);length=int(r['protein_length']) if mode=='protein' else int(r['region_length'])
  flags=dict(support_dual_workflow=dual,support_close_at_least_3=int(r['close_num'])>=3,support_bacterial_archaeal_best_distant=ba,support_length_at_least_150=length>=150,support_entropy_at_most_2_5=float(r['entropy'])<=2.5);extra=sum(flags.values())
  ctier=tier(int(r['contig_length']),int(r['genome_N50']),int(r['dist_to_edge']),float(r['hgt_contig_fraction']),extra)
  composition=r.get('contig_non_eukaryotic_fraction','');clean=composition!='' and float(composition)<=.5
  gate=yes(r['pass_mapping']) if mode=='protein' else yes(r['pass_manuscript_coverage_variant']) and int(r['internal_stop_codons'])<2 and int(r['protein_length'])>=80
  eligible=gate and clean and ctier.startswith('Tier');r.update(**flags,extra_support_n=extra,contig_context_tier=ctier,eligible_for_phylogeny=eligible,dual_workflow_supporting_loci=';'.join(x['locus_id'] for x in evidence),dual_workflow_overlap_definition='reciprocal_at_least_50pct_genomic_overlap_same_strand',tier3_support_version=2)
  if mode=='protein':r.update(pass_paper_contig_tier=ctier.startswith('Tier'),high_confidence_flank_context=eligible and yes(r['bilateral_eukaryotic_flanks']))
  if old!=eligible:changes.append(dict(workflow=r.get('_mode',mode),locus_id=r['locus_id'],genome=r['genome'],accession=q,old_eligible=old,new_eligible=eligible,outside_old_catalog_locus=r['outside_old_catalog_locus'],**flags))
for mode,context,mapping in configs:
 rr=[dict((k,v) for k,v in r.items() if not k.startswith('_')) for r in protein if r['_mode']==mode];md=O/'discovery'/mapping
 for name,selected in [('all_protein_candidate_loci.tsv',rr),('protein_loci_for_phylogeny.tsv',[r for r in rr if yes(r['eligible_for_phylogeny'])])]:
  p=md/name
  if p.exists() and not p.with_suffix('.before_full_support.tsv').exists():shutil.copy2(p,p.with_suffix('.before_full_support.tsv'))
  write(p,selected)
 bystem=defaultdict(list)
 for r in rr:bystem[Path(r['flank_fasta']).parent].append(r)
 for folder,selected in bystem.items():write(folder/'loci.tsv',selected)
 summary=json.loads((md/'protein_mapping.complete.json').read_text());summary.update(eligible_loci=sum(yes(r['eligible_for_phylogeny']) for r in rr),novel_eligible_loci=sum(yes(r['eligible_for_phylogeny']) and yes(r['outside_old_catalog_locus']) for r in rr),tier3_support_version=2);(md/'protein_mapping.complete.json').write_text(json.dumps(summary,indent=2))
nd=O/'discovery/nucleotide_context';write(nd/'all_expanded_region_metrics.tsv',nucleotide);eligible=[r for r in nucleotide if yes(r['eligible_for_phylogeny'])];write(nd/'nucleotide_loci_for_phylogeny.tsv',eligible)
s=json.loads((nd/'nucleotide_context.complete.json').read_text());s.update(eligible_loci=len(eligible),novel_eligible_loci=sum(yes(r['outside_old_catalog_locus']) for r in eligible),tier3_support_version=2);(nd/'nucleotide_context.complete.json').write_text(json.dumps(s,indent=2))
write(D/'tier3_full_support_changes.tsv',changes);(D/'all_contig_support.complete.json').write_text(json.dumps({'five_supports_evaluated':True,'changed_loci':len(changes),'newly_eligible_loci':sum(r['new_eligible'] for r in changes)},indent=2));print('Full support pass',len(changes),'changed loci',flush=True)
# Existing completed trees are reused; only newly eligible query trees are added.
for script in ['run_protein_phylogeny.py','audit_finished_protein_trees.py','run_recovered_protein_phylogeny.py','run_missing_protein_phylogeny.py']:
 with (O/'logs'/('final_update_'+script+'.log')).open('w') as log:subprocess.run(['python','-u',str(O/'scripts'/script)],stdout=log,stderr=subprocess.STDOUT,check=True)
missing=json.loads((O/'discovery/missing_protein_context/protein_context.complete.json').read_text());(O/'discovery/missing_protein_context/followup.complete.json').write_text(json.dumps({'proteins_checked':missing['queried_proteins'],'distribution_pass':missing['passes_distribution'],'genomic_context_and_phylogeny_complete':True}))
(D/'protein_tree_updates.complete.json').write_text(json.dumps({'all_protein_routes_updated_after_full_contig_support':True}));print('Additional protein tree evaluation complete',flush=True)
