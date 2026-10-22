from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,re,json,hashlib
import numpy as np
csv.field_size_limit(10000000)
ROOT=Path(__file__).resolve().parents[5];R=ROOT/'analysis/revision_pnas/revision_round2_assessment';O=Path(__file__).resolve().parents[1];T=O/'tables';G=R/'innovation_execution/genome_integration';F=R/'innovation_final8989_20260915';A=ROOT/'align/protein_align/genome_divide'
manifest={}
def rows(p):
 p=Path(p);manifest[str(p)]=dict(bytes=p.stat().st_size)
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as h:yield from csv.DictReader(h,delimiter='\t')
def write(name,rr,fields=None):
 rr=list(rr)
 with (T/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields or (list(rr[0]) if rr else ['none']),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def toks(x):return {z for z in re.split('[;,|]',x or '') if z not in {'','NA','nan','-'}}
def yes(x):return str(x).lower() in {'true','1'}
panel={r['species']:r for r in rows(R/'supplementary_tables_reader_20261006/supporting_data/S2/Q_species_counts.tsv') if r['scenario']=='augmented_same_reference'}
assert len(panel)==48
cat={r['hgt_id']:{k:r[k] for k in ['hgt_id','Specie','Genome','Accession','Data_type','final_addition']} for r in rows(F/'inputs/final_catalogue_8989.tsv')}
assert len(cat)==8989
ledger={r['record_id']:r for r in rows(G/'candidate_record_ledger.tsv') if r['record_id'] in cat}
inv={r['unit_id']:r for r in rows(R/'manuscript_evidence_completion_20261003/data/final48_locus_inventory.tsv.gz')}
assert sum(r['role']=='candidate' for r in inv.values())==2299
candset={u for u,r in inv.items() if r['role']=='candidate'}
# Verify the stored plotting counts against the actual model inventory.
countchecks=[]
for sp,r in panel.items():
 for role,fc,nc in [('candidate','C','N'),('background','B','M')]:
  z=[a for a in inv.values() if a['species']==sp and a['role']==role];cc=len({a['OG'] for a in z});nn=len(z)
  countchecks.append(dict(species=sp,role=role,expected_families=int(float(r[fc])),observed_families=cc,expected_loci=int(float(r[nc])),observed_loci=nn,passed=cc==int(float(r[fc])) and nn==int(float(r[nc]))))
assert all(r['passed'] for r in countchecks);write('01_count_reconstruction.tsv',countchecks)
# Within-species CD-HIT output, with observed representative/member identities.
clusters={};membermap=defaultdict(list);filechecks=[]
def readclusters(p):
 cluster=None;members=[]
 with p.open() as f:
  for line in f:
   if line.startswith('>Cluster'):
    if cluster is not None:yield cluster,members
    cluster=line.strip()[1:];members=[]
   else:
    m=re.search(r'\d+\s+(\d+)(aa|nt), >(.+)\.\.\. (.+)',line.strip())
    if not m:raise ValueError((p,line))
    ident=100. if m[4]=='*' else float(re.search(r'([\d.]+)%',m[4])[1]);members.append(dict(accession=m[3],length=int(m[1]),identity_pct=ident,is_rep=m[4]=='*'))
 if cluster is not None:yield cluster,members
for p in sorted(A.glob('GCA_*.HGT_verify_filtered.entropy.clustered.fa.clstr')):
 genome=p.name.split('.HGT_verify')[0];manifest[str(p)]=dict(bytes=p.stat().st_size);n=0;nm=0;ids=[]
 for c,mm in readclusters(p):
  reps=[m for m in mm if m['is_rep']];assert len(reps)==1
  rep=reps[0]['accession'];key=(genome,rep);assert key not in clusters;clusters[key]=mm;n+=1;nm+=len(mm)
  for m in mm:membermap[genome,m['accession']].append(rep)
  ids.extend(m['identity_pct'] for m in mm if not m['is_rep'])
 filechecks.append(dict(genome=genome,clusters=n,member_records=nm,nonrepresentatives=nm-n,min_observed_member_identity=min(ids) if ids else '',source=str(p)))
write('02_per_species90_files.tsv',filechecks)
# Trace final records into both saved 50%-named outputs; do not assume either is the catalogue filter.
globalrows=[];maps50={}
for name in ['105_parasite_prot_HGT.entropy_clu50_cluster.tsv','clu50_cluster.tsv']:
 p=A/name;manifest[str(p)]=dict(bytes=p.stat().st_size);pairs=list(csv.reader(p.open(),delimiter='\t'));m2r=defaultdict(set)
 for rep,mem in pairs:m2r[mem].add(rep)
 maps50[name]=m2r
 for r in cat.values():
  if ':' in r['Accession']:continue
  reps=m2r.get(r['Accession'],set());globalrows.append(dict(map_file=name,record_id=r['hgt_id'],species=r['Specie'],accession=r['Accession'],in_map=bool(reps),is_representative=r['Accession'] in reps,representatives=';'.join(sorted(reps)),in_48_species=r['Specie'] in panel))
write('03_final_records_vs50.tsv',globalrows)
# Accessions of all members of clusters anchored to a final catalogue record in the 48-species comparison.
edges=[];needed=defaultdict(set)
for rid,r in cat.items():
 sp=r['Specie'];ac=r['Accession'];genome=r['Genome']
 if sp not in panel or ':' in ac:continue
 for rep in membermap.get((genome,ac),[]):
  for m in clusters[genome,rep]:
   needed[sp].add(m['accession']);edges.append(dict(anchor_record=rid,species=sp,clade=panel[sp]['clade'],genome=genome,anchor_accession=ac,representative=rep,member_accession=m['accession'],member_identity_pct=m['identity_pct'],member_is_representative=m['is_rep'],anchor_is_representative=ac==rep))
for r in cat.values():
 if r['Specie'] in panel and ':' not in r['Accession']:needed[r['Specie']].add(r['Accession'])
print('prepared',len(inv),'inventory units',len(edges),'cluster edges',sum(map(len,needed.values())),'needed accession keys',flush=True)
og={r['unit_id']:r['inferred_family'] for r in rows(G/'refinement_20260914/unit_orthogroup_membership.tsv.gz')}
# Preserve every different model; use original mapping fields, never transfer HGT status via a cluster.
accunit=defaultdict(set);units={};ridunit=defaultdict(set);scan=[]
keys=['species','unit_id','gene_id','sequence_id','contig','strand','coding_start_0','coding_end_0','coding_structure','eligible_structure','analysis_role','aa_length_measured']
for inp in rows(G/'species_inputs.tsv'):
 sp=inp['species']
 if sp not in panel:continue
 p=G/'measurements'/(inp['species_id']+'.tsv.gz')
 if not p.exists():continue
 cnt=0;matched=0
 for a in rows(p):
  cnt+=1;u=a['unit_id'];acs=toks(a['accessions'])|{a['accession']};hit=acs&needed[sp]
  if hit or u in candset:
   z={k:a.get(k,'') for k in keys};z['OG']=og.get(u,'');z['current_role']=inv.get(u,{}).get('role','not_in_panel');z['source_type']='original_measurement';z['accessions']=';'.join(sorted(acs));units[u]=z
   for ac in hit:accunit[sp,ac].add(u)
   if u in candset:
    for rid in toks(a['direct_catalogue_records'])|toks(a['unique_candidate_records']):
     if rid in cat:ridunit[rid].add(u)
   matched+=1
 scan.append(dict(species=sp,measured_models_scanned=cnt,retained_for_audit=matched,source=str(p)))
# Added nematode reference reconstruction; exactly the input used for final counts.
for a in rows(R/'metazoa_extension_20260922/tables/gene_units_recovered.tsv.gz'):
 sp=a['species'];u=a['unit_id']
 if sp not in panel or sp in {x['species'] for x in scan}:continue
 z={k:a.get(k,'') for k in keys};z['OG']=a['OG'];z['current_role']=inv.get(u,{}).get('role','not_in_panel');z['source_type']='metazoan_recovery';z['accessions']='';units[u]=z
 if u in candset:
  for rid in toks(a['record_ids']):
   if rid in cat:ridunit[rid].add(u)
# Recovered models lack accession lists; join through the existing original sequence aliases.
newsp=set(panel)-{x['species'] for x in scan};seq_units=defaultdict(set)
for u,a in units.items():
 if a['species'] in newsp:seq_units[a['species'],a['sequence_id']].add(u)
if newsp:
 for a in rows(R/'metazoa_extension_20260922/tables/original_protein_input_annotations.tsv.gz'):
  sp=a['species'];ac=a['accession']
  if sp in newsp and ac in needed[sp]:accunit[sp,ac].update(seq_units.get((sp,a['sequence_id']),set()))
# Final strict ledger is a second independent linkage source, restricted to actual candidate units.
for r in rows(F/'inputs/final8989_record_locus_family_ledger.tsv'):
 for u in toks(r['strict_own_unit_ids'])&candset:ridunit[r['record_id']].add(u)
write('04_scan_coverage.tsv',scan);write('05_relevant_model_inventory.tsv',units.values())
linked={u for v in ridunit.values() for u in v};print('candidate units with record link',len(linked),flush=True)
# Cross-table evidence, including cases not included in the actual 48-species locus panel.
annot=[];potential={};statuses=Counter()
for e in edges:
 sp=e['species'];ac=e['member_accession'];au=ridunit.get(e['anchor_record'],set());mu=accunit.get((sp,ac),set());mu= {u for u in mu if units[u]['species']==sp};same=set(mu)&au
 sameog={u for u in mu if any(units[u]['OG'] and units[u]['OG']==inv[a]['OG'] for a in au)}
 good={u for u in sameog if yes(units[u]['eligible_structure']) and u not in au}
 distinct={u for u in good if all((units[u]['contig'],units[u]['strand'],units[u]['coding_structure'])!=(inv[a]['contig'],inv[a]['strand'],inv[a]['coding_structure']) for a in au)}
 # Distinct coding structures may overlap; only non-overlapping gene spans form the strict scenario.
 nonoverlap={u for u in distinct if all(units[u]['contig']!=inv[a]['contig'] or int(units[u]['coding_end_0'])<=int(inv[a]['coding_start_0']) or int(units[u]['coding_start_0'])>=int(inv[a]['coding_end_0']) for a in au)}
 uncounted={u for u in nonoverlap if u not in inv}
 status='anchor_not_in_strict_panel' if not au else 'no_model_link' if not mu else 'already_counted_candidate' if mu&candset else 'same_gene_as_anchor' if same else 'ambiguous_multiple_member_models' if len(mu)>1 else 'potential_uncounted_distinct_locus' if uncounted else 'other_model_or_family_or_eligibility'
 # Only single model per member and a single anchor locus. No automatic acceptance.
 eligible=uncounted if len(mu)==1 and len(au)==1 else set()
 for u in eligible:
  potential.setdefault(u,dict(units[u],anchor_records=set(),member_accessions=set(),cluster_representatives=set()))
  potential[u]['anchor_records'].add(e['anchor_record']);potential[u]['member_accessions'].add(ac);potential[u]['cluster_representatives'].add(e['representative'])
 annot.append(dict(**e,anchor_in_strict_panel=bool(au),anchor_units=';'.join(sorted(au)),member_model_count=len(mu),member_units=';'.join(sorted(mu)),member_roles=';'.join(sorted({units[u]['current_role'] for u in mu})),same_family_eligible_distinct_units=';'.join(sorted(nonoverlap)),sensitivity_extra_units=';'.join(sorted(eligible)),status=status));statuses[status]+=1
write('06_cluster90_member_to_locus.tsv',annot)
for u,z in potential.items():
 for k in ['anchor_records','member_accessions','cluster_representatives']:z[k]=';'.join(sorted(z[k]))
 z['HGT_confirmed_by_this_audit']=False
write('07_potential_uncounted_loci_NOT_NEW_HGT.tsv',potential.values())
# Model sequence redundancy is audited in both actual roles, independent of candidate clustering.
blocks=defaultdict(list)
for u,a in inv.items():blocks[a['species'],a['role'],a['sequence_id']].append(a)
identical=[]
for (sp,role,sid),rr in blocks.items():
 if len(rr)>1:identical.append(dict(species=sp,role=role,sequence_id=sid,models=len(rr),distinct_coding_structures=len({(x['contig'],x['strand'],x['coding_structure']) for x in rr}),units=';'.join(x['unit_id'] for x in rr)))
write('08_identical_sequence_distinct_models.tsv',identical)
# Frozen n=3 draws reproduce the published descriptive contrasts; adding unconfirmed loci is a sensitivity only.
spcounts=[]
for sp,r in panel.items():
 extra=[z for z in potential.values() if z['species']==sp];spcounts.append(dict(species=sp,clade=r['clade'],candidate_families=int(float(r['C'])),background_families=int(float(r['B'])),candidate_loci=int(float(r['N'])),background_loci=int(float(r['M'])),potential_additional_loci=len(extra),conditional_candidate_loci=int(float(r['N']))+len(extra)))
write('09_species_counts_and_sensitivity.tsv',spcounts)
cts=[]
for scen in ['published','all_eligible90_members_assumed_HGT_SENSITIVITY_ONLY']:
 rng=np.random.default_rng(20260922);arrays={}
 for cl in ['Amoebozoa','Apicomplexa','Euglenozoa','Metazoa']:
  z=sorted([r for r in spcounts if r['clade']==cl],key=lambda r:r['species']);arr=[]
  for i in range(1000):
   zz=[z[k] for k in rng.choice(len(z),3,replace=False)];c=sum(v['candidate_families'] for v in zz);b=sum(v['background_families'] for v in zz);n=sum(v['candidate_loci'] if scen=='published' else v['conditional_candidate_loci'] for v in zz);m=sum(v['background_loci'] for v in zz);arr.append((c/b,(n/c)/(m/b),n/m))
  arrays[cl]=np.asarray(arr)
 for cl in ['Amoebozoa','Apicomplexa','Euglenozoa']:
  ratio=arrays[cl]/arrays['Metazoa']
  for j,metric in enumerate(['family_ratio','within_family_loci_ratio','locus_ratio']):
   vals=ratio[:,j];cts.append(dict(scenario=scen,clade=cl,metric=metric,fold=np.exp(np.log(vals).mean()),q025=np.quantile(vals,.025),q975=np.quantile(vals,.975),formal_P_value=False))
write('10_descriptive_sensitivity.tsv',cts)
# Map 50% clusters onto OGs only when current locus links are reliable; neither clustering nor OG implies an HGT event.
mapmain=maps50['105_parasite_prot_HGT.entropy_clu50_cluster.tsv'];fifty=[]
for rid,uu in ridunit.items():
 r=cat[rid]
 for u in uu:
  for rep in mapmain.get(r['Accession'],set()):fifty.append(dict(record_id=rid,species=r['Specie'],unit_id=u,OG=inv[u]['OG'],accession=r['Accession'],cluster50_rep=rep,is_rep=r['Accession']==rep))
write('11_strict_panel_50cluster_OG_links.tsv',fifty)
summary=dict(final_records=len(cat),panel_species=len(panel),candidate_loci=len(candset),background_loci=sum(a['role']=='background' for a in inv.values()),candidate_units_with_record_link=len(linked),cluster90_files=len(filechecks),cluster90_members=sum(r['member_records'] for r in filechecks),cluster90_reps=sum(r['clusters'] for r in filechecks),edge_status_counts=dict(statuses),potential_additional_loci=len(potential),potential_by_clade=dict(Counter(panel[z['species']]['clade'] for z in potential.values())),global50_maps={},candidate_identical_sequence_blocks=sum(r['role']=='candidate' for r in identical),background_identical_sequence_blocks=sum(r['role']=='background' for r in identical))
for name in maps50:
 rr=[r for r in globalrows if r['map_file']==name];summary['global50_maps'][name]=dict(final_protein_records=len(rr),mapped=sum(r['in_map'] for r in rr),mapped_nonrepresentatives=sum(r['in_map'] and not r['is_representative'] for r in rr))
(O/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(O/'input_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
