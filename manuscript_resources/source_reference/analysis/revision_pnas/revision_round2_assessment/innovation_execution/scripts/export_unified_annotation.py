"""Project unified annotations to accession/species; retain no-hit and unresolved locus states."""
from pathlib import Path
import pandas as pd,json,ast,re,gzip,csv,collections
O=Path(__file__).resolve().parents[1];A=O/'annotation';R=O.parent/'innovation_function_ecology';V=O.parents[1]/'virus_supplement'
status=json.loads((A/'run_status.json').read_text());assert status.get('complete'), 'Annotation not complete: no partial export'
f=Path(status['annotation_file']);header=json.loads((A/'annotation_columns.json').read_text());ann={}
with f.open() as h:
 for line in h:
  if line.startswith('#query'):header=line.lstrip('#').rstrip().split('\t')
  elif not line.startswith('#') and line.strip():
   assert header;vals=line.rstrip('\n').split('\t');assert len(vals)==len(header);r=dict(zip(header,vals));assert r['query'] not in ann;ann[r['query']]=r
# Freeze and reuse original descriptive module definitions without executing the old inventory builder.
source=V/'scripts/build_innovation_inventory.py';tree=ast.parse(source.read_text());rules=None
for n in tree.body:
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='modules' for t in n.targets):rules=ast.literal_eval(n.value)
assert rules;compiled={k:re.compile(v,re.I) for k,v in rules.items()};(A/'uniform_module_rules.json').write_text(json.dumps(rules,indent=2))
fields=['Description','Preferred_name','KEGG_ko','KEGG_Pathway','PFAMs','eggNOG_OGs','EC','GOs','CAZy']
for r in ann.values():
 text=' '.join(r.get(k,'') for k in fields);r['functional_modules']=';'.join(k for k,p in compiled.items() if p.search(text))
 r['has_informative_annotation']=any(r.get(k,'') not in ['','-','NA'] for k in ['Description','KEGG_ko','PFAMs'])
protein={}
with (V/'inputs/protein_to_genome.tsv').open() as h:
 for r in csv.DictReader(h,delimiter='\t'):
  if r['sequence_available'].lower()=='true':assert r['accession'] not in protein;protein[r['accession']]=r
# Match exact genome+protein ID only. Same sequence is not proof of same locus.
genes=collections.defaultdict(set)
with (O/'GFF_all_gene_models.tsv').open() as h:
 for r in csv.DictReader(h,delimiter='\t'):
  if r['protein_id'] and r['gene_id']:genes[(r['genome'],r['protein_id'])].add(r['gene_id'])
candidates=pd.read_csv(R/'combined_candidate_function_inventory.tsv',sep='\t',keep_default_na=False)
needed=set(zip(candidates.genome,candidates.accession));candidate_maps={};summary=collections.defaultdict(collections.Counter);modules=collections.defaultdict(collections.Counter);gene_modules=collections.defaultdict(set);n=0
out=A/'all105_protein_annotations.tsv.gz';temp=Path(str(out)+'.partial')
keys=['accession','sequence_id','species','genome','aa_length','sha256','annotation_status','gene_mapping_status','gene_id','functional_modules','has_informative_annotation']+[k for k in header if k!='query']
with gzip.open(temp,'wt',newline='') as h,(A/'sequence_aliases.tsv').open() as al:
 w=csv.DictWriter(h,fieldnames=keys,delimiter='\t');w.writeheader()
 for s in csv.DictReader(al,delimiter='\t'):
  p=protein[s['accession']];a=ann.get(s['sequence_id']);g=genes.get((p['genome'],s['accession']),set());gm='unique_exact_GFF_gene' if len(g)==1 else 'multiple_GFF_genes' if g else 'unresolved_alias_or_missing_GFF'
  r={k:'' for k in keys};r.update(s,species=p['species'],genome=p['genome'],annotation_status='annotated' if a else 'no_transferred_annotation',gene_mapping_status=gm,gene_id=';'.join(sorted(g)))
  if a:r.update({k:v for k,v in a.items() if k in r})
  w.writerow(r);n+=1;summary[p['species']]['input_accessions']+=1;summary[p['species']][r['annotation_status']]+=1;summary[p['species']][gm]+=1
  if a and a['has_informative_annotation']:summary[p['species']]['informative_annotation']+=1
  for mod in (a['functional_modules'].split(';') if a else []):
   if not mod:continue
   modules[(p['species'],mod)]['annotated_accessions']+=1
   if len(g)==1:gene_modules[(p['species'],mod)].add(next(iter(g)))
  if (p['genome'],s['accession']) in needed:candidate_maps[(p['genome'],s['accession'])]=r
assert n==len(protein)==1352152
assert len(summary)==105;temp.replace(out)
pd.DataFrame([dict(species=s,**v) for s,v in sorted(summary.items())]).fillna(0).to_csv(A/'species_annotation_coverage.tsv',sep='\t',index=False)
pd.DataFrame([dict(species=s,module=m,annotated_accessions=modules[(s,m)]['annotated_accessions'],resolved_GFF_gene_units=len(gene_modules[(s,m)]),counts_are_not_complete_biological_copies=True) for s in sorted(summary) for m in rules]).to_csv(A/'species_uniform_module_background.tsv',sep='\t',index=False)
rows=[]
for _,c in candidates.iterrows():
 a=candidate_maps.get((c.genome,c.accession));r={k:c[k] for k in ['record_id','scope','locus_unit','species','genome','accession','data_type','tree_class','latest_tree_outside_whitelist']}
 r.update(uniform_mapping_status='exact_proteome_accession' if a else 'unresolved_full_CDS_or_alias',uniform_annotation_status=a['annotation_status'] if a else 'not_evaluated',uniform_modules=a['functional_modules'] if a else '',gene_mapping_status=a['gene_mapping_status'] if a else '',gene_id=a['gene_id'] if a else '',sequence_id=a['sequence_id'] if a else '')
 for k in fields:r['uniform_'+k]=a.get(k,'') if a else ''
 rows.append(r)
cd=pd.DataFrame(rows);cd.to_csv(A/'candidate_uniform_annotation_crosswalk.tsv',sep='\t',index=False)
report={'complete':True,'input_accessions':n,'unique_sequence_annotations':len(ann),'species':len(summary),'candidate_exact_proteome_accessions':int(cd.uniform_mapping_status.eq('exact_proteome_accession').sum()),'candidate_unresolved_full_CDS_or_alias':int(cd.uniform_mapping_status.ne('exact_proteome_accession').sum()),'module_rules':'historical_v1_regex_descriptive_not_curated_orthogroups','gene_unit_counts':'resolved GFF subset only, not total copy number','no_annotation_is_not_function_absence':True,'export_file':str(out)}
(A/'uniform_export_status.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

import subprocess
subprocess.run(["python",str(O/"ecology_hgt_analysis/build_current_report.py")],check=True)
