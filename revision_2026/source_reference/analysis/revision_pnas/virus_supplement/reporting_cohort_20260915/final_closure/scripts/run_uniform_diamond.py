"""New uniform DIAMOND review authorized by user; do not censor weak Close hits.
Only the fixed priority19 query records are searched here, not a de novo 105-species screen.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from Bio.SeqIO.FastaIO import SimpleFastaParser
import csv,json,subprocess,hashlib,time
F=Path(__file__).resolve().parents[1];D=F.parent;B=D.parent;S=F/'diamond_review';S.mkdir(exist_ok=True)
rows=list(csv.DictReader((D/'tables/viral_priority19_all_evidence.tsv').open(),delimiter='\t'));ids={r['locus_id'] for r in rows}
with (B/'downstream_novel/queries_all.faa').open() as h:records=[(t.split()[0],s) for t,s in SimpleFastaParser(h) if t.split()[0] in ids]
assert len(records)==19
query=S/'priority19.faa';query.write_text(''.join('>'+t+'\n'+s+'\n' for t,s in records));qhash=hashlib.sha256(query.read_bytes()).hexdigest()
plan=dict(scope='fixed_19_priority_records_only',not_de_novo_105_species_screen=True,query_sha256=qhash,search='DIAMOND blastp sensitive; evalue1e-10; max-target-seqs2000; no minimum identity or coverage at search stage',grouping='original recipient/host hierarchy; virus classification restored; ambiguous multi-source taxids remain flagged',initial_gate='original Out proportion/count/AI/first-Out identity and coverage; species entropy threshold held fixed from full historical distribution, not estimated from19',entropy_calibration='fixed legacy whole-species viral-repaired cutoff; explicitly a frozen external calibration for new DIAMOND review',sampling='original length>=80aa identity>=20%; original dynamic quotas and extraction; preserve low-coverage Close hits',output_scope='nr_only and nr_plus_deduplicated_viral_panel; not literal reproduction of original MMseqs',saturated_query_action='flag for expanded hit-cap recheck before certification',original8985='frozen; no deletion')
(S/'analysis_plan.json').write_text(json.dumps(plan,indent=2))
(S/'software_version.txt').write_text(subprocess.check_output(['diamond','version'],text=True))
def run(panel,db,threads):
 out=S/(panel+'.tsv');status=S/(panel+'.status.json');cmd=['diamond','blastp','--query',str(query),'--db',str(db),'--out',str(out)+'.part','--threads',str(threads),'--sensitive','--evalue','1e-10','--max-target-seqs','2000','--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']
 fp=dict(query_sha256=qhash,command=cmd,database_bytes=Path(str(db)+'.dmnd' if not str(db).endswith('.dmnd') else db).stat().st_size)
 if status.exists() and out.exists():
  prev=json.load(status.open())
  if prev.get('status')=='complete' and prev.get('fingerprint')==fp:return prev
 status.write_text(json.dumps(dict(status='running',fingerprint=fp),indent=2));start=time.time()
 with (S/(panel+'.log')).open('w') as h:res=subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT)
 z=dict(status='complete' if res.returncode==0 else 'failed',fingerprint=fp,elapsed_seconds=time.time()-start,returncode=res.returncode)
 if res.returncode==0:Path(str(out)+'.part').replace(out)
 status.write_text(json.dumps(z,indent=2));print(panel,z['status'],z['elapsed_seconds'],flush=True);return z
with ThreadPoolExecutor(max_workers=2) as pool:
 jobs=[pool.submit(run,'nr',Path('EXTERNAL_DATA_ROOT/Data/DB/nr/diamond_nr/nr_rep_seq.dmnd'),24),pool.submit(run,'viral',B/'databases/combined_viral.dmnd',2)]
 results=[f.result() for f in as_completed(jobs)]
(S/'search_complete.json').write_text(json.dumps(dict(complete=all(r['status']=='complete' for r in results),jobs=results),indent=2))
