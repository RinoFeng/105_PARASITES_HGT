"""Run minimally patched original classifiers followed by unmodified original filters."""
from pathlib import Path
import subprocess,csv,json,concurrent.futures,time,sys,os
D=Path(__file__).resolve().parents[1];R=D.parents[3]
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
species=list(csv.reader((R/'105_name2accfil_withTaxid.txt').open(),delimiter='\t'))
configs=[('protein','align/protein_align/genome_divide','align/protein_align/genome_divide/target_nt_taxid_lineage.txt','HGT_verify_virus_fixed.prot.py','HGT_filter_entropy.prot.py'),('nucleotide_r1','align/nucl_align/genome_divided','align/nucl_align/genome_divided/target_nt_taxid_lineage.txt','HGT_verify_virus_fixed.nucl.py','HGT_filter_entropy.py'),('nucleotide_r2','align/second_align','align/target_nt_taxid_lineage2.txt','HGT_verify_virus_fixed.nucl.py','HGT_filter_entropy.py')]
def one(task):
 (sp,g,tid),config=task;mode,folder,tax,ver,fil=config;outdir=D/'classification_only'/mode;outdir.mkdir(parents=True,exist_ok=True);marker=outdir/(g+'.status.json');raw=outdir/(g+'.fixed_metrics.tsv');passed=outdir/(g+'.fixed_pass.tsv');source=R/folder/(g+'.acc2taxid.out');log=outdir/(g+'.log')
 if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':return json.loads(marker.read_text())
 cmd1=[sys.executable,str(D/'scripts'/ver),str(R/tax),str(source),tid];cmd2=[sys.executable,str(D/'scripts'/fil),str(raw),'-o',str(passed)]
 status=dict(species=sp,genome=g,workflow=mode,status='running',source=str(source),source_bytes=source.stat().st_size,commands=[cmd1,cmd2],start=time.strftime('%Y-%m-%d %H:%M:%S'));marker.write_text(json.dumps(status,indent=2))
 with log.open('w') as lf,raw.with_suffix('.part').open('w') as ff:cp=subprocess.run(cmd1,stdout=ff,stderr=lf)
 if cp.returncode==0:
  raw.with_suffix('.part').replace(raw)
  with log.open('a') as lf:cp=subprocess.run(cmd2,stdout=lf,stderr=subprocess.STDOUT)
 status.update(returncode=cp.returncode,end=time.strftime('%Y-%m-%d %H:%M:%S'),status='complete' if cp.returncode==0 and passed.exists() else 'failed')
 if status['status']=='complete':
  with passed.open() as f:status['fixed_pass_count']=sum(1 for _ in f)-1
 marker.write_text(json.dumps(status,indent=2));print(mode,sp,status['status'],status.get('fixed_pass_count',''),flush=True);return status
if __name__=='__main__':
 tasks=[(s,c) for c in configs for s in species];results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
  for r in pool.map(one,tasks):results.append(r)
 (D/'reports/classification_only_status.json').write_text(json.dumps(results,indent=2));print('FINISHED',len(results),'failed',sum(r['status']!='complete' for r in results),flush=True)
