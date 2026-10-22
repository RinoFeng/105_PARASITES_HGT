"""Contig composition from original window summaries overlaid with viral-fixed replay.
Unassigned windows remain unassigned. Coordinates are 1-based inclusive.
"""
import csv,json,time,concurrent.futures
from pathlib import Path
from collections import defaultdict,Counter
O=Path(__file__).resolve().parents[1];R=O.parents[2];D=O/'discovery/contig_evidence';D.mkdir(exist_ok=True)
def coord(q):
 if ':' in q:
  c,x=q.rsplit(':',1);a,b=x.split('-')
 else:c,a,b=q.rsplit('-',2)
 return c,int(a),int(b)
def vote(x):
 out=float(x['Out_pct']);close=int(x['Close_num']);distant=int(x.get('Distant_num',0));ai=float(x['alien_index']);selfn=int(x['Self_num']);lineage=x.get('First_distant_taxinfo','')
 close_pct=close/(close+distant) if close+distant else 0
 non_euk=ai>0 and out>.8 and any(n in lineage.split(';') for n in ['Bacteria','Archaea','Viruses'])
 euk=ai<0 or close_pct>.2 or selfn>1
 # Preserve both flags: conflicting assignments are not evidence of clean recipient context.
 return 'conflicting' if non_euk and euk else 'non_eukaryotic' if non_euk else 'eukaryotic' if euk else 'unassigned'
def one(s):
 g=s['genome'];fixed=O/'recovery/rounded_legacy/nucleotide_r1'/(g+'.viral_and_changed.tsv');marker=fixed.with_name(g+'.summary.json')
 while not marker.exists():time.sleep(15)
 replacements={r['query']:r for r in csv.DictReader(fixed.open(),delimiter='\t')}
 old=R/'align/nucl_align/genome_divided'/(g+'.HGT_verify.out');counts=defaultdict(Counter);seen=set();rows=[]
 dest=D/(g+'.window_votes.tsv')
 with dest.open('w') as out:
  w=csv.writer(out,delimiter='\t');w.writerow(['genome','contig','start_1based','end_1based','query','vote','viral_fix_overlay'])
  if old.exists():
   with old.open() as f:
    next(f)
    for r in csv.DictReader(f,delimiter='\t'):
     q=r['Accession'];seen.add(q);x=replacements.get(q,r);c,a,b=coord(q);v=vote(x);counts[c][v]+=1;w.writerow([g,c,a,b,q,v,q in replacements])
  for q,x in replacements.items():
   if q in seen:continue
   c,a,b=coord(q);v=vote(x);counts[c][v]+=1;w.writerow([g,c,a,b,q,v,True])
 for c,ct in counts.items():
  assign=ct['non_eukaryotic']+ct['eukaryotic'];conservative_non=ct['non_eukaryotic']+ct['conflicting'];conservative_all=assign+ct['conflicting']
  rows.append({'genome':g,'contig':c,'eukaryotic_windows':ct['eukaryotic'],'non_eukaryotic_windows':ct['non_eukaryotic'],'conflicting_windows':ct['conflicting'],'unassigned_windows':ct['unassigned'],'assignable_windows_excluding_conflicts':assign,'non_eukaryotic_fraction_excluding_conflicts':ct['non_eukaryotic']/assign if assign else '', 'non_eukaryotic_fraction_counting_conflicts':conservative_non/conservative_all if conservative_all else '', 'contamination_flag_conservative':conservative_non/conservative_all>.5 if conservative_all else '', 'evidence_basis':'stored_original_window_summary_with_viral_fixed_overlay'})
 print(s['species'],len(rows),flush=True);return rows
if __name__=='__main__':
 specs=list(csv.DictReader((O/'inputs/species_manifest.tsv').open(),delimiter='\t'));rows=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for result in pool.map(one,specs):rows.extend(result)
 with (D/'all105_contig_votes.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 (D/'contig_evidence.complete.json').write_text(json.dumps({'genomes':len(specs),'contigs_with_old_alignment_summary':len(rows),'conflicting_window_policy':'explicit flags; conservative fraction counts conflicts as non-eukaryotic; unassigned excluded'}))
