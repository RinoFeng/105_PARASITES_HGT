"""Prepare locus-resolved comparative units. Non-catalogue never means proven native."""
from pathlib import Path
import csv,gzip,json,collections,hashlib
import pandas as pd
from Bio import SeqIO
D=Path(__file__).resolve().parent;O=D.parent;A=O/'annotation';R=O.parent/'innovation_function_ecology'
c=pd.read_csv(R/'combined_candidate_function_inventory.tsv',sep='\t',keep_default_na=False)
links=collections.defaultdict(set)
for _,r in c.iterrows():
 for a in [r.accession]+r.cluster_member_accessions.split(','):
  if a and a not in ['NA','nan']:links[(r.genome,a.strip())].add(r.record_id)
# Independent full-CDS/GFF audit can identify aliases and coordinate candidates. Keep evidence type.
p=pd.read_csv(O/'priority_fullCDS_GFF_crosswalk.tsv',sep='\t',keep_default_na=False)
priority=collections.defaultdict(set)
for _,r in p.iterrows():
 if r.mapping_status=='unique_GFF_match' and r.GFF_gene_id:
  priority[(r.genome,r.GFF_gene_id)].add(r.record_id)
units={};candidate_links=[];excluded=collections.Counter()
keys=['accession','sequence_id','species','genome','aa_length','gene_id','functional_modules','COG_category','Description','Preferred_name','KEGG_ko','PFAMs','EC','eggNOG_OGs']
with gzip.open(A/'all105_protein_annotations.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['gene_mapping_status']!='unique_exact_GFF_gene':excluded[r['gene_mapping_status']]+=1;continue
  key=(r['genome'],r['gene_id']);ids=set(links.get((r['genome'],r['accession']),set()));alias_ids=priority.get(key,set())
  for cid in ids|alias_ids:candidate_links.append((cid,r['species'],r['genome'],r['gene_id'],r['accession'],'catalogue_accession_or_cluster_member' if cid in ids else 'priority_unique_GFF_match'))
  if key not in units:
   units[key]={k:r[k] for k in keys};units[key]['accessions']=set();units[key]['catalogue_records']=set();units[key]['direct_catalogue_records']=set()
  u=units[key];u['accessions'].add(r['accession']);u['catalogue_records'].update(ids|alias_ids);u['direct_catalogue_records'].update(ids)
  if (int(r['aa_length']),r['accession'])>(int(u['aa_length']),u['accession']):
   for k in keys:u[k]=r[k]
# All copies at distinct annotated loci retained even when sequence-identical.
rows=[];seq_to_units=collections.defaultdict(list)
for n,key in enumerate(sorted(units),1):
 u=units[key];u['unit_id']=f'GENE{n:07d}';u['catalogue_status']='catalogue_linked' if u['catalogue_records'] else 'non_catalogue_origin_unverified'
 for k in ['accessions','catalogue_records','direct_catalogue_records']:u[k]=';'.join(sorted(u[k]))
 rows.append(u);seq_to_units[u['sequence_id']].append(u['unit_id'])
pd.DataFrame(rows).to_csv(D/'resolved_gene_units.tsv.gz',sep='\t',index=False,compression='gzip')
pd.DataFrame(candidate_links,columns=['record_id','species','genome','gene_id','accession','mapping_basis']).drop_duplicates().to_csv(D/'catalogue_gene_links.tsv',sep='\t',index=False)
n=0
with (D/'gene_representatives.faa').open('w') as out:
 for s in SeqIO.parse(A/'unique_sequences.faa','fasta'):
  for uid in seq_to_units.get(s.id,[]):out.write('>'+uid+'\n'+str(s.seq)+'\n');n+=1
assert n==len(units)
report={'gene_units':n,'species_with_resolved_units':len({u['species'] for u in rows}),'catalogue_linked_units':sum(u['catalogue_status']=='catalogue_linked' for u in rows),'mapped_catalogue_records':len({r[0] for r in candidate_links}),'excluded_accessions':dict(excluded),'isoform_rule':'longest protein at exact same GFF gene; deterministic accession tie-break','non_catalogue_is_not_native':True,'sequence_identity_does_not_merge_distinct_loci':True,'biological_copy_number_fully_validated':False}
(D/'gene_units_status.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
