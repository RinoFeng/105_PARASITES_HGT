"""Close canonical protein gates, retaining alternate ORF records separately.
The old catalogue is frozen. Missing evidence never counts as a pass.
"""
from pathlib import Path
from collections import defaultdict,Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import csv,json,subprocess,shutil,hashlib
F=Path(__file__).resolve().parents[1];D=F.parent;B=D.parent;R=B.parents[2];A=B/'harmonized_final_20260915';S=F/'diamond_review'
for n in ['dedup90','dedup80']:(F/n).mkdir(exist_ok=True)
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def fa(p):
 with p.open() as h:return [(t.split()[0],s) for t,s in SimpleFastaParser(h)]
def write(p,rs):
 if not rs:return
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(dict.fromkeys(k for r in rs for k in r)),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def fasta(p,rs):
 with p.open('w') as h:
  for k,s in rs:h.write('>'+k+'\n'+s+'\n')
records=rows(D/'tables/viral_priority19_all_evidence.tsv');byid={r['locus_id']:r for r in records};pre=rows(F/'tables/diamond_pre_tree_gate19.tsv');source_review={r['locus_id']:r for r in rows(F/'tables/viral_source_topology_review.tsv')};query=dict(fa(B/'downstream_novel/queries_all.faa'));canonical={r['locus_id']:r for r in rows(F/'tables/record_to_canonical_locus.tsv')};categories=json.load((B/'reports/original_tree_category_map.json').open())['categories'];tree={};products={}
for arm in ['trees_diamond','trees_diamond_C20']:
 p=F/arm/'tree_classification/tree_classification.tsv';tree[arm]={r['hgt_id']:r for r in rows(p)} if p.exists() else {}
 for p in (F/arm/'jobs').glob('*/status.json'):products[arm,p.parent.name]=json.load(p.open())
# Early pool is declared, frozen historical initial pass pool plus newly passing
# canonical protein queries in this supplementary review (before downstream gates).
initial=defaultdict(list)
for p in pre:
 if p['canonical_role']=='canonical_protein_location' and p['initial_distribution']=='PASS':initial[p['genome']].append(byid[p['locus_id']])
poolids={};required=set()
for g,rr in initial.items():
 poolids[g]={r['Accession'] for r in rows(A/'route_filter_validation/protein'/(g+'.fixed_pass.tsv'))}|{r['accession'] for r in rr};required.update(poolids[g])
seq=[];seen={};alias=[];duplicate=[];source=Path('EXTERNAL_DATA_ROOT/HGT/parasite/new_157_parasite/protein/sequence.fasta')
if required:
 with source.open() as h:
  for t,s in SimpleFastaParser(h):
   token=t.split()[0];k=token if token in required else (token.split('|')[1] if token.startswith(('sp|','tr|')) and len(token.split('|'))>2 else token)
   if k not in required:continue
   if k in seen:duplicate.append(k);continue
   seen[k]=s;seq.append((k,s))
   if k!=token:alias.append(dict(accession=k,original_token=token,sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
# New DIAMOND query FASTA is an authorized source for supplement-only entries.
for g,rr in initial.items():
 for r in rr:
  if r['accession'] not in seen:seen[r['accession']]=query[r['locus_id']];seq.append((r['accession'],query[r['locus_id']]))
early={};earlymanifest=[]
for g,rr in initial.items():
 ids=poolids[g];missing=ids-set(seen);bad=set(duplicate)&ids;inp=F/'dedup90'/(g+'.fa');fasta(inp,[(k,s) for k,s in seq if k in ids]);out=F/'dedup90'/(g+'.clustered.fa');cmd=['cd-hit','-c','0.9','-i',str(inp),'-o',str(out)]
 with (F/'logs'/(g+'.final90.log')).open('w') as h:z=subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT)
 groups=[];group=[];rep={}
 if z.returncode==0:
  for line in Path(str(out)+'.clstr').read_text().splitlines():
   if line.startswith('>'):
    if group:groups.append(group)
    group=[]
   else:group.append((line.split('>',1)[1].split('...',1)[0],line.endswith('*')))
  if group:groups.append(group)
  for group in groups:
   rk=next(k for k,star in group if star)
   for k,_ in group:rep[k]=rk
 earlymanifest.append(dict(genome=g,expected_n=len(ids),extracted_n=len(ids)-len(missing),missing_accessions=';'.join(sorted(missing)),duplicate_accessions=';'.join(sorted(bad)),command=json.dumps(cmd),returncode=z.returncode,input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),scope='frozen_historical_repaired_initial_pool_plus_this_supplement_initial_pass_queries;not_a_full_DIAMOND_reclassification_of_105_proteomes'))
 for r in rr:
  k=r['accession'];eq=seen.get(k,'').rstrip('*')==query[r['locus_id']].rstrip('*');rk=rep.get(k,'');early[r['locus_id']]=dict(representative=rk,status='PENDING_incomplete_or_conflicting_pool' if missing or bad or not eq or not rk else 'PASS_complete_declared_pool_representative' if rk==k else 'FAIL_early90_nonrepresentative',input_sequence_matches=eq)
write(F/'tables/final_early90_pool_manifest.tsv',earlymanifest);write(F/'tables/final_early90_aliases.tsv',alias)
base_stages=['initial_distribution','source_priority','BUSCO','animal_conservation','EGT_symbiont','genome_presence','original_mapping_contig','revision_contig_tier','bilateral_flanks','reported_region_sequence_consistency','hit_cap','source_ambiguity']
stages=base_stages+['early90','tree_sampling','tree_primary','tree_C20_sensitivity','viral_source_consistency','final80','old_catalog_locus_overlap']
post=[]
for p in pre:
 z=dict(p);q=p['locus_id'];er=early.get(q,{});z['early90']=er.get('status','NOT_REACHED');z['early90_representative']=er.get('representative','');z['tree_sampling']='NOT_REACHED';z['tree_primary']='NOT_REACHED';z['tree_C20_sensitivity']='NOT_REACHED'
 for arm,label in [('trees_diamond','tree_primary'),('trees_diamond_C20','tree_C20_sensitivity')]:
  t=tree[arm].get(q);prod=products.get((arm,q),{})
  if t:
   valid=prod.get('products_valid',False) and prod.get('ufboot_convergence') is not None and prod['ufboot_convergence']>=.99
   z[label]='PASS_original_whitelist_and_QC' if valid and t['primary_class'] in categories and t['long_branch']=='False' else 'FAIL_tree_category_or_long_branch_or_convergence' if valid else 'PENDING_tree_products_or_convergence'
   z[label+'_class']=t['primary_class'];z[label+'_flags']=t['flags'];z[label+'_long_branch']=t['long_branch'];z[label+'_convergence']=prod.get('ufboot_convergence');z[label+'_file']=str(F/arm/'jobs'/q/'iqtree.treefile')
  if label=='tree_primary' and t:
   prep={r['hgt_id']:r for r in rows(F/'tables/diamond_tree_input_status.tsv')}.get(q,{})
   z['tree_sampling']='PASS_original_sampler_extractor_with_identical_query_legacy_context' if prep.get('status')=='ok' and not prep.get('missing_subjects') else 'PENDING_missing_or_conflicting_references'
 z['viral_source_consistency']=source_review.get(q,{}).get('viral_source_consistency','NOT_REACHED');z['final80']='NOT_REACHED';z['old_catalog_locus_overlap']='NOT_REACHED';z['admissible_before_final_dedup']=z['canonical_role']=='canonical_protein_location' and all(z[k].startswith('PASS') for k in stages if k not in ['final80','old_catalog_locus_overlap']);post.append(z)
# Compare against every available old coordinate and preserve missing-coordinate coverage.
old=rows(A/'inputs/historical8985_frozen.tsv');oldkeys={(r['Genome'],r['Accession']) for r in old};positions=defaultdict(list);missing_old=defaultdict(int)
for r in old:
 try:a,b=sorted([int(float(r['hgt_start'])),int(float(r['hgt_end']))]);c=r['Contig_Name'];assert c and a>0;positions[r['Genome'],c].append((a,b,r['hgt_id']))
 except (KeyError,ValueError,AssertionError):missing_old[r['Genome']]+=1
accepted_groups=defaultdict(list)
for z in post:
 if z['admissible_before_final_dedup']:accepted_groups[z['genome']].append(z)
cl=R/'analysis/Results_article_new/Remove_EGT/rm_wrong_contig/Analysis_Clustering_prot/step1_intra_species';manifest80=[]
for g,rr in accepted_groups.items():
 species=byid[rr[0]['locus_id']]['species'];src=cl/(species.replace(' ','_').replace('/','_')+'.in.fasta');data=fa(src);have=dict(data);conflict=set()
 for z in rr:
  r=byid[z['locus_id']];k=r['accession'];s=query[r['locus_id']]
  if k in have and have[k].rstrip('*')!=s.rstrip('*'):conflict.add(k)
  if k not in have:data.append((k,s));have[k]=s
 inp=F/'dedup80'/(g+'.fa');fasta(inp,data);prefix=F/'dedup80'/g;tmp=F/'dedup80'/(g+'.tmp');cmd=['mmseqs','easy-cluster',str(inp),str(prefix),str(tmp),'--min-seq-id','0.8','-c','0.8','--cov-mode','1','--threads','4','-v','1']
 with (F/'logs'/(g+'.final80.log')).open('w') as h:res=subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT)
 rep={};members=defaultdict(list)
 if res.returncode==0:
  for rk,k in csv.reader(Path(str(prefix)+'_cluster.tsv').open(),delimiter='\t'):rep[k]=rk;members[rk].append(k)
 manifest80.append(dict(genome=g,source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),input_n=len(data),returncode=res.returncode,command=json.dumps(cmd),sequence_conflicts=';'.join(sorted(conflict)),old_coordinate_missing_n=missing_old[g]))
 for z in rr:
  r=byid[z['locus_id']];k=r['accession'];rk=rep.get(k,'');mem=members.get(rk,[]);oldmem=[a for a in mem if (g,a) in oldkeys]
  z['final80_representative']=rk;z['final80_members']=';'.join(mem);z['final80_old8985_members']=';'.join(oldmem)
  z['final80']='PENDING_missing_or_conflicting_pool' if not rk or k in conflict or len(data)!=len(have) else 'FAIL_already_represented_in_old8985' if oldmem else 'PASS_species_representative' if rk==k else 'FAIL_final80_nonrepresentative'
  a,b=int(r['start_1based']),int(r['end_1based']);over=[hid for c,d,hid in positions[g,r['contig']] if max(0,min(b,d)-max(a,c)+1)>=.5*min(b-a+1,d-c+1)];z['overlapping_old_ids']=';'.join(over)
  z['old_catalog_locus_overlap']='FAIL_overlaps_old_catalog_locus' if over else 'PASS_no_overlap_old_catalog_coordinates' if missing_old[g]==0 else 'PENDING_old_coordinate_coverage_incomplete'
for z in post:
 q=z['locus_id'];z['all_required_stages_passed']=z['canonical_role']=='canonical_protein_location' and all(z[k].startswith('PASS') for k in stages);z['add_to_primary']=z['all_required_stages_passed'];z['nonpass_stages']=';'.join(k for k in stages if not z[k].startswith('PASS'))
 if z['canonical_role']!='canonical_protein_location':z['final_status']='supporting_alternate_ORF_not_independent_addition'
 elif z['all_required_stages_passed']:z['final_status']='certified_supplementary_candidate'
 elif any(z[k].startswith('FAIL') for k in stages):z['final_status']='not_added_failed_required_gate_or_represented'
 else:z['final_status']='unresolved_not_added'
write(F/'tables/final_candidate_gate19.tsv',post);write(F/'tables/final80_pool_manifest.tsv',manifest80)
accepted=[z['locus_id'] for z in post if z['add_to_primary']];unresolved=[z['locus_id'] for z in post if z['final_status']=='unresolved_not_added'];summary=dict(scope='fixed19_records_17_canonical_protein_locations',original8985_unchanged=True,accepted_locus_ids=accepted,certified_additions=len(accepted),certified_catalogue_size=8985+len(accepted),unresolved_locus_ids=unresolved,final_count_determined_within_declared_scope=not unresolved,status_counts=dict(Counter(z['final_status'] for z in post)),calibration='frozen_external_historical_whole_species_entropy_threshold',not_a_full_105_proteome_reanalysis=True,nt_alternate_records_not_independently_certified=True)
(F/'final_gate_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
