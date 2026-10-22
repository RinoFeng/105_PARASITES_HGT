"""Map viral-associated proteins to their own assemblies; retain every reported locus.
Contig tiers follow the submitted SI. Flank composition is an additional evidence field.
"""
import csv,json,time,subprocess,concurrent.futures
from pathlib import Path
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
O=Path(__file__).resolve().parents[1];R=O.parents[2];D=O/'discovery/protein_genome_mapping';D.mkdir(exist_ok=True)
def wait(p):
 while not p.exists():time.sleep(15)
def table(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def tier(length,n50,edge,fraction,extra):
 if fraction>.7:return 'excluded_region_fraction'
 if length>=100000 and edge>=10000:return 'Tier1_strict_100kb_middle'
 if length>=max(20000,min(50000,n50)) and edge>=5000:return 'Tier2_adaptive'
 if length>=20000 and edge>=2000 and extra>=3:return 'Tier3_extra_support'
 return 'excluded_contig_context'
def one(s):
 g=s['genome'];rr=BY.get(g,[])
 if not rr:return []
 sd=D/Path(s['window_fasta']).stem;sd.mkdir(exist_ok=True);qf=sd/'queries.faa'
 with qf.open('w') as f:
  for r in rr:f.write('>'+r['accession']+'\n'+SEQ[r['accession']]+'\n')
 gff=sd/'mapped.gff'
 if not (sd/'mapping.done.json').exists():
  cmd=['miniprot','-I','-t','2','--gff','--outc=0.8',s['genome_path'],str(qf)]
  with (sd/'mapped.gff.part').open('w') as f,(sd/'miniprot.log').open('w') as log:subprocess.run(cmd,stdout=f,stderr=log,check=True)
  (sd/'mapped.gff.part').replace(gff);(sd/'mapping.done.json').write_text(json.dumps({'command':cmd}))
 lookup={r['accession']:r for r in rr};mrnas=[];exons=defaultdict(int)
 for line in gff.open():
  if line.startswith('#'):continue
  z=line.rstrip('\n').split('\t')
  if len(z)!=9:continue
  a=dict(x.split('=',1) for x in z[8].split(';') if '=' in x)
  if z[2]=='mRNA':mrnas.append((z,a))
  elif z[2]=='CDS':exons[a.get('Parent','')]+=1
 votes=defaultdict(list)
 for x in table(O/'discovery/contig_evidence'/(g+'.window_votes.tsv')):votes[x['contig']].append((int(x['start_1based']),int(x['end_1based']),x['vote']))
 seqs={t.split()[0]:seq for t,seq in SimpleFastaParser(open(s['genome_path']))};rows=[];mapped=set()
 with (sd/'candidate_regions_with_10kb_flanks.fna').open('w') as flanks:
  for i,(z,a) in enumerate(mrnas,1):
   target=a['Target'].split();q=target[0];r=lookup[q];mapped.add(q);c=z[0];start=int(z[3]);end=int(z[4]);length=len(seqs[c]);edge=min(start-1,length-end);frac=(end-start+1)/length
   cov=(int(target[2])-int(target[1])+1)/len(SEQ[q]);identity=float(a.get('Identity',0));frameshift=int(a.get('Frameshift',0))
   extra=int(int(r['close_num'])>=3)+int(int(r['protein_length'])>=150)+int(float(r['entropy'])<=2.5)
   ctier=tier(length,int(s['genome_N50']),edge,frac,extra);cv=CV.get((g,c),{});fraction=cv.get('non_eukaryotic_fraction_counting_conflicts','')
   left=sum(v=='eukaryotic' and a0>=max(1,start-10000) and b0<start for a0,b0,v in votes[c]);right=sum(v=='eukaryotic' and a0>end and b0<=min(length,end+10000) for a0,b0,v in votes[c])
   intersect=[h['hgt_id'] for h in FINAL.get((g,c),[]) if max(0,min(end,int(h['hgt_end']))-max(start,int(h['hgt_start']))+1)>=.5*min(end-start+1,int(h['hgt_end'])-int(h['hgt_start'])+1)]
   locus=f'VP_{sd.name}_{i:06d}';fa=max(1,start-10000);fb=min(length,end+10000);flanks.write(f'>{locus}|{c}:{fa}-{fb}|candidate={start}-{end}|strand={z[6]}\n{seqs[c][fa-1:fb]}\n')
   coordok=cov>=.9 and identity>=.9;clean=fraction!='' and float(fraction)<=.5
   rows.append(dict(r,locus_id=locus,contig=c,start_1based=start,end_1based=end,strand=z[6],mapping_identity=identity,mapping_query_coverage=cov,mapping_frameshifts=frameshift,mapping_cds_segments=exons[a['ID']],contig_length=length,genome_N50=s['genome_N50'],dist_to_edge=edge,hgt_contig_fraction=frac,contig_context_tier=ctier,extra_support_n=extra,contig_non_eukaryotic_fraction=fraction,left_eukaryotic_flank_windows=left,right_eukaryotic_flank_windows=right,overlapping_old_hgt_ids=';'.join(intersect),outside_old_catalog_locus=not intersect and r['old_catalog_member']=='False',pass_mapping=coordok,pass_paper_contig_tier=ctier.startswith('Tier'),pass_contig_composition=clean,bilateral_eukaryotic_flanks=left>0 and right>0,eligible_for_phylogeny=coordok and ctier.startswith('Tier') and clean,high_confidence_flank_context=coordok and ctier.startswith('Tier') and clean and left>0 and right>0,flank_fasta=str(sd/'candidate_regions_with_10kb_flanks.fna'),inference_status='candidate_requires_phylogeny'))
 with (sd/'unmapped_accessions.txt').open('w') as f:
  for q in sorted(set(lookup)-mapped):f.write(q+'\n')
 if rows:
  with (sd/'loci.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 print(s['species'],'candidates',len(rr),'loci',len(rows),'tree eligible',sum(x['eligible_for_phylogeny'] for x in rows),flush=True);return rows
if __name__=='__main__':
 for p in [O/'discovery/protein_context/protein_context.complete.json',O/'discovery/contig_evidence/contig_evidence.complete.json']:wait(p)
 specs=table(O/'inputs/genome_window_manifest.tsv');BY=defaultdict(list)
 for r in table(O/'discovery/protein_context/protein_candidates_for_genome_context.tsv'):BY[r['genome']].append(r)
 need={r['accession'] for rr in BY.values() for r in rr};SEQ={}
 for p in [O/'inputs/all_105_proteins.faa',O/'inputs/supplemental_recovered_proteins.faa',O/'inputs/late_recovered_alias_proteins.faa']:
  with p.open() as f:
   for t,seq in SimpleFastaParser(f):
    if t.split()[0] in need:SEQ[t.split()[0]]=seq
 CV={(r['genome'],r['contig']):r for r in table(O/'discovery/contig_evidence/all105_contig_votes.tsv')};FINAL=defaultdict(list)
 for r in table(R/'analysis/revision_pnas/HGT_candidate_verify/main_phylogeny_retained_HGT_8985.tsv'):FINAL[(r['Genome'],r['Contig_Name'])].append(r)
 results=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
  for rows in pool.map(one,specs):results.extend(rows)
 for name,rows in [('all_protein_candidate_loci.tsv',results),('protein_loci_for_phylogeny.tsv',[r for r in results if r['eligible_for_phylogeny']])]:
  with (D/name).open('w') as f:
   if results:
    w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 (D/'protein_mapping.complete.json').write_text(json.dumps({'candidate_queries':len(need),'mapped_loci':len(results),'eligible_loci':sum(r['eligible_for_phylogeny'] for r in results),'novel_eligible_loci':sum(r['eligible_for_phylogeny'] and r['outside_old_catalog_locus'] for r in results)},indent=2))
