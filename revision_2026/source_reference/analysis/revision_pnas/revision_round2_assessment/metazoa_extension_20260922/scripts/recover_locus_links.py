from pathlib import Path
from collections import defaultdict,Counter
import sys,re,json,gzip,hashlib,warnings
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
B=Path(__file__).resolve().parents[1];R=B.parent;O=R/'innovation_execution';G=O/'genome_integration';F=R/'innovation_final8989_20260915';T=B/'tables';sys.path.insert(0,str(G))
from measure import parse_gff,model_metrics
inputs={}
def read(p,**kw):
 p=Path(p);inputs[str(p.resolve())]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,sep='\t',keep_default_na=False,low_memory=False,**kw)
def write(n,r):pd.DataFrame(r).to_csv(T/n,sep='\t',index=False,compression='infer')
def chunks_for_species(p,cols,ss):
 p=Path(p);inputs[str(p.resolve())]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.concat([t[t.species.isin(ss)] for t in pd.read_csv(p,sep='\t',usecols=cols,keep_default_na=False,chunksize=150000,low_memory=False)],ignore_index=True)
def normset(x):return {z.strip() for z in re.split('[;,|]',x) if z.strip() not in {'','-','NA','nan'}}
refs=read(T/'reference_manifest.tsv');ss=set(refs.species);cols=['accession','sequence_id','species','genome','aa_length','sha256','COG_category','KEGG_ko','PFAMs','eggNOG_OGs']
ann=chunks_for_species(O/'annotation/all105_protein_annotations.tsv.gz',cols,ss);write('original_protein_input_annotations.tsv.gz',ann)
fi=chunks_for_species(O/'ecology_macro_execution/inferred_families/family_input_membership.tsv.gz',['species','sequence_id','inferred_family'],ss);ogmap=fi.groupby(['species','sequence_id']).inferred_family.agg(lambda z:set(z)).to_dict()
final=read(F/'inputs/final_catalogue_8989.tsv',usecols=['hgt_id']);final=set(final.hgt_id);historical=read(G/'candidate_record_ledger.tsv');historical=historical[historical.species.isin(ss)]
blocked=set(read(F/'macro/excluded_ordinary_background_OGs.tsv').OG);global_added_blocked=set();units=[];modelrows=[];audit=[];summ=[];bridges=[]
for ref in refs.itertuples():
 sp=ref.species;aa=ann[ann.species==sp];hr=historical[historical.species==sp];byhash=defaultdict(list);byacc=defaultdict(set);byseq={}
 for x in aa.to_dict('records'):
  byhash[x['sha256']].append(x);byacc[x['accession']].add(x['sequence_id']);byseq.setdefault(x['sequence_id'],x)
 directseq=defaultdict(set);possibleseq=defaultdict(set);coord=[]
 for x in hr.itertuples():
  for sid in byacc.get(x.accession,set()):directseq[sid].add(x.record_id);possibleseq[sid].add(x.record_id)
  for ac in normset(x.cluster_member_accessions):
   for sid in byacc.get(ac,set()):possibleseq[sid].add(x.record_id)
  ma=re.fullmatch(r'([^:]+):(\d+)-(\d+)',x.accession)
  if ma:coord.append((ma[1],min(int(ma[2]),int(ma[3]))-1,max(int(ma[2]),int(ma[3])),x.record_id))
 # Public protein headers and full-AA equality are used only as version-link evidence.
 with gzip.open(ref.protein_path,'rt') as h:
  for p in SeqIO.parse(h,'fasta'):
   seq=str(p.seq).upper();sha=hashlib.sha256(seq.encode()).hexdigest();matches=byhash.get(sha,[])
   rr={rid for z in matches for rid in directseq.get(z['sequence_id'],set())}
   for rid in sorted(rr):bridges.append(dict(record_id=rid,species=sp,reference_assembly=ref.reference_assembly,published_protein_id=p.id,full_AA_identical=True,identity_pct=100,query_coverage=1,reference_coverage=1,same_original_reference=ref.same_original_reference,locus_validated=False))
 features,genes,unresolved=parse_gff(ref.gff_path);op=gzip.open if str(ref.reference_genome_path).endswith('.gz') else open
 with op(ref.reference_genome_path,'rt') as h:genome={r.id:str(r.seq).upper() for r in SeqIO.parse(h,'fasta')}
 allmodels=defaultdict(list);possiblegenes=defaultdict(set);directgenes=defaultdict(set);eligible_by_record=defaultdict(set);selected_by_gene={};allcoords=defaultdict(list)
 for co,st,en,rid in coord:allcoords[co].append((st,en,rid))
 stats=Counter()
 for gid,choices in genes.items():
  gf=features[gid];touch=set()
  if ref.same_original_reference:
   for st,en,rid in allcoords.get(gf['contig'],[]):
    if gf['start']<en and gf['end']>st:touch.add(rid);possiblegenes[gid].add(rid)
  for cid,parts,ancestors,flags in choices:
   stats['CDS_models']+=1;co=parts[0]['contig'];strand=parts[0]['strand'];ordered=sorted(parts,key=lambda p:p['start'],reverse=strand=='-')
   if co not in genome:stats['contig_not_found']+=1;continue
   if ordered[0]['phase'] not in {'0','1','2'}:stats['initial_phase_unknown']+=1;continue
   fr=[genome[co][p['start']:p['end']] for p in ordered]
   if strand=='-':fr=[str(Seq(t).reverse_complement()) for t in fr]
   cds=''.join(fr)[int(ordered[0]['phase']):];tables={int(p['attrs']['transl_table']) for p in parts+ancestors if p['attrs'].get('transl_table')};table=next(iter(tables)) if len(tables)==1 else 1
   try:protein=str(Seq(cds).translate(table=table,cds=True))
   except Exception:
    if len(cds)%3:stats['nontriplet_no_complete_translation']+=1;continue
    with warnings.catch_warnings():warnings.simplefilter('ignore');protein=str(Seq(cds).translate(table=table))
    if protein.endswith('*'):protein=protein[:-1]
   if '*' in protein:stats['internal_stop']+=1;continue
   sha=hashlib.sha256(protein.upper().encode()).hexdigest();hits=byhash.get(sha,[])
   if not hits:stats['not_exact_original_protein']+=1;continue
   sids={p['sequence_id'] for p in hits}
   if len(sids)!=1:stats['ambiguous_original_sequence_identity']+=1;continue
   sid=next(iter(sids));mm=model_metrics(parts,ancestors,genome[co],sha)
   if flags and 'annotation_flags' in mm:mm['annotation_flags']=';'.join(sorted((set(mm['annotation_flags'].split(';'))-{''})|set(flags)))
   caseids=directseq.get(sid,set());possibleids=possibleseq.get(sid,set())|touch;possiblegenes[gid]|=possibleids
   for rid in caseids:directgenes[rid].add(gid)
   if mm.get('model_error'):stats['invalid_model']+=1;continue
   geneok=gf['contig']==co and gf['strand']==strand and 0<=gf['start']<gf['end']<=len(genome[co]) and mm['coding_start_0']>=gf['start'] and mm['coding_end_0']<=gf['end']
   eligible=geneok and mm['translation_exact'] and mm['complete_start_stop'] and mm['phase_consistent'] and mm['initial_phase']==0 and mm['cds_ambiguous']==0 and not mm['annotation_flags'] and not mm['short_CDS_gap_under20']
   fams=ogmap.get((sp,sid),set());fam=next(iter(fams)) if len(fams)==1 else ''
   row=dict(species=sp,reference_assembly=ref.reference_assembly,same_original_reference=ref.same_original_reference,gene_id=gid,CDS_ID=cid,contig=co,strand=strand,sequence_id=sid,OG=fam,OG_mapping_status='unique_original_OG' if len(fams)==1 else 'missing_or_ambiguous_original_OG',aa_sha256=sha,own_records=';'.join(sorted(caseids)),possible_records=';'.join(sorted(possibleids)),eligible_structure=bool(eligible),gene_bounds_valid=bool(geneok),**mm)
   allmodels[gid].append(row);modelrows.append(row);stats['exact_models']+=1
  choices=allmodels.get(gid,[])
  if not choices:continue
  # Same representative rule as the original recovery; exact coding ambiguity remains excluded.
  best=max(choices,key=lambda m:(m['complete_start_stop'],m['aa_length_measured'],m['CDS_ID']));sig={(z['contig'],z['strand'],z['coding_structure']) for z in choices if z['sequence_id']==best['sequence_id']}
  best=dict(best);best['unique_representative_structure']=len(sig)==1;best['eligible_structure']=best['eligible_structure'] and len(sig)==1;selected_by_gene[gid]=best
  for rid in normset(best['own_records']):
   if best['eligible_structure']:eligible_by_record[rid].add(gid)
 # Match original own-record strict selection: both overall direct locus and eligible locus unique.
 strict_record={};speciesaudit=[]
 for x in hr.itertuples():
  ds=directgenes.get(x.record_id,set());es=eligible_by_record.get(x.record_id,set());strictgid=next(iter(es)) if len(ds)==len(es)==1 else None
  if strictgid:state='own_sequence_exact_eligible_representative'
  elif len(ds)>1:state='multiple_possible_gene_loci'
  elif ds:state='own_sequence_match_not_eligible_representative'
  elif any(x.record_id in ids for ids in possiblegenes.values()):state='coordinate_or_alias_only_not_strict'
  elif x.data_type=='Nucleotide' and not ref.same_original_reference:state='old_coordinates_on_different_reference_not_transferred'
  else:state='no_full_sequence_locus_link'
  if strictgid and x.record_id in final:strict_record[x.record_id]=strictgid
  speciesaudit.append(dict(record_id=x.record_id,species=sp,in_final_catalogue=x.record_id in final,original_accession=x.accession,original_data_type=x.data_type,original_tree_class=x.tree_class,original_contig_tier=x.contig_tier,reference_assembly=ref.reference_assembly,same_original_reference=ref.same_original_reference,original_contig_evidence_transferred=bool(ref.same_original_reference),direct_genes=';'.join(sorted(ds)),eligible_representative_genes=';'.join(sorted(es)),strict_gene_id=strictgid or '',status=state,new_HGT_event=False))
 audit+=speciesaudit;bygid=defaultdict(set)
 for rid,gid in strict_record.items():bygid[gid].add(rid)
 added=set()
 for gid,rr in selected_by_gene.items():
  if possiblegenes.get(gid) and rr['OG']:added.add(rr['OG'])
 global_added_blocked |= added
 localblock=blocked|added;out=[]
 for gid,rr in selected_by_gene.items():
  seq=rr['sequence_id'];an=byseq[seq];fam=rr['OG'];rid=bygid.get(gid,set());role='candidate' if rid and fam else 'ordinary_background' if rr['eligible_structure'] and fam and not possiblegenes.get(gid) and fam not in localblock else 'excluded'
  uid='GENE_'+hashlib.sha256((ref.original_genome+'|'+gid if ref.same_original_reference else ref.reference_assembly+'|'+gid).encode()).hexdigest()[:14]
  row=dict(rr,unit_id=uid,role=role,record_ids=';'.join(sorted(rid)),historical_possible_records=';'.join(sorted(possiblegenes.get(gid,set()))),COG_category=an['COG_category'],eggnog_present=an['eggNOG_OGs'] not in {'','-','NA','nan'},KO_present=an['KEGG_ko'] not in {'','-','NA','nan'},PFAM_present=an['PFAMs'] not in {'','-','NA','nan'},original_annotation_exactly_reused=True)
  units.append(row);out.append(row)
 stats.update(dict(exact_genes=len(out),eligible_genes=sum(z['eligible_structure'] for z in out),strict_final_records=len(strict_record),strict_candidate_loci=sum(z['role']=='candidate' for z in out),strict_candidate_OGs=len({z['OG'] for z in out if z['role']=='candidate'}),ordinary_background_loci=sum(z['role']=='ordinary_background' for z in out),ordinary_background_OGs=len({z['OG'] for z in out if z['role']=='ordinary_background'})))
 summ.append(dict(species=sp,reference_assembly=ref.reference_assembly,same_original_reference=ref.same_original_reference,original_protein_observations=aa.sequence_id.nunique(),GFF_genes_with_CDS=len(genes),unresolved_GFF_parent_models=unresolved,final_records=sum(x['in_final_catalogue'] for x in speciesaudit),**stats));print(sp,dict(stats),flush=True)
 del genome
# Conservatively quarantine any newly recovered historical candidate family across all recovered species.
allblocked=blocked|global_added_blocked
for z in units:
 if z['role']=='ordinary_background' and z['OG'] in allblocked:z['role']='excluded'
write('gene_units_recovered.tsv.gz',units);write('all_exact_CDS_models.tsv.gz',modelrows);write('candidate_locus_audit.tsv',audit);write('protein_version_bridges.tsv',bridges);write('species_recovery_summary.tsv',summ);write('new_historical_blocked_OGs.tsv',pd.DataFrame({'OG':sorted(global_added_blocked-blocked)}));write('all_historical_blocked_OGs.tsv',[dict(OG=o) for o in sorted(allblocked)]);write('input_manifest.tsv',[dict(path=p,sha256=h) for p,h in inputs.items()]);
status=dict(status='complete',species=len(refs),records_audited=len(audit),strict_records=sum(z['strict_final_records'] for z in summ),newly_quarantined_OGs=len(global_added_blocked-blocked),source_labels_changed=False,tree_ratings_changed=False,new_HGT_records=0)
(B/'RECOVERY_STATUS.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status))
