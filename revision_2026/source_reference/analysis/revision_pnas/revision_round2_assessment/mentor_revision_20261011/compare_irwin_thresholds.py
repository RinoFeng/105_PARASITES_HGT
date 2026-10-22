from pathlib import Path
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
import csv, io, json, zipfile, hashlib, subprocess

ROOT = Path('PROJECT_ROOT')
B = ROOT/'analysis/revision_pnas/virus_supplement'
O = Path(__file__).resolve().parent/'irwin_thresholds'
O.mkdir(exist_ok=True)
def read(p):
    with p.open() as f: return list(csv.DictReader(f, delimiter='\t'))
def write(p, rows):
    with p.open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta(p):
    with p.open() as f: return {t.split()[0]:s.upper().rstrip('*') for t,s in SimpleFastaParser(f)}
lit=read(B/'reports/irwin_direction_resolved_20260915/S2_1423_direction_resolved.tsv')
lm={x['irwin_table2_id']:x for x in lit}
seq={}
with zipfile.ZipFile(B/'literature/irwin2022/Irwin_et_al_VirusHGT.zip') as z:
    with z.open('dryad/1.Datasets_and_annotations/final_clustering_euk_virus_prok.fasta') as f:
        for t,s in SimpleFastaParser(io.TextIOWrapper(f)):
            sid=t.split()[0]
            if sid in lm: seq[sid]=s.upper().rstrip('*')
assert len(seq)==1423
(O/'author_1423.faa').write_text(''.join('>'+s+'\n'+seq[s]+'\n' for s in sorted(seq)))
cat=read(B/'reporting_cohort_20260915/final_closure/tables/final_catalogue_8989.tsv')
oldaa=fasta(B/'inputs/final_6878_protein_records.faa');oldnt=fasta(B/'inputs/final_2107_nucleotide_records.fna')
newaa=fasta(B/'reporting_cohort_20260915/final_closure/accepted_new.proteins.faa')
manifest=read(B/'literature/irwin2022/comparison/query_manifest.tsv')
gt={x['genome']:x['species_taxid'] for x in manifest}
parents={};ranks={}
with (B/'databases/taxonomy/nodes.dmp').open() as f:
    for line in f:
        a=line.split('|');parents[a[0].strip()]=a[1].strip();ranks[a[0].strip()]=a[2].strip()
def species(t):
    seen=set()
    while t not in seen:
        if ranks.get(t)=='species':return t
        seen.add(t);t=parents.get(t,'1')
    return t
qm={};qseq={};aa=[];nt=[]
for i,r in enumerate(cat):
    q='Q%05d'%i
    typ='nt' if r['Data_type']=='Nucleotide' else 'aa'
    s=oldnt[r['Accession']] if typ=='nt' else newaa[r['hgt_id']] if r['hgt_id'] in newaa else oldaa[r['Accession']]
    qseq[q]=s;qm[q]=dict(query=q,record_id=r['hgt_id'],accession=r['Accession'],species=r['Specie'],species_taxid=gt[r['Genome']],representation=typ)
    (nt if typ=='nt' else aa).append('>'+q+'\n'+s+'\n')
(O/'candidates_aa.faa').write_text(''.join(aa));(O/'candidates_nt.fna').write_text(''.join(nt));write(O/'query_manifest.tsv',list(qm.values()))
commands=[]
def run(cmd):
    commands.append(cmd)
    p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    if p.returncode: raise RuntimeError(p.stderr[-3000:])
run(['diamond','makedb','--in',str(O/'author_1423.faa'),'-d',str(O/'author_1423')])
for mode,fn in [('blastp','candidates_aa.faa'),('blastx','candidates_nt.fna')]:
    run(['diamond',mode,'--query',str(O/fn),'--db',str(O/'author_1423.dmnd'),'--out',str(O/(mode+'.tsv')),'--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','--evalue','1e-10','--id','30','--max-target-seqs','0','--sensitive','--threads','8'])
    print(mode,'finished',flush=True)
(O/'commands.json').write_text(json.dumps(commands,indent=2))
hits=[]
for mode in ['blastp','blastx']:
    with (O/(mode+'.tsv')).open() as f:
        for z in csv.reader(f,delimiter='\t'):
            q,s=z[:2];m=qm[q];l=lm[s]
            hits.append(dict(**m,irwin_id=s,irwin_species=l['study_species'],contamination_score=float(l['contamination_score']),same_species=m['species_taxid']==species(s.split('.')[0]),identity_pct=float(z[2]),query_coverage=(abs(int(z[5])-int(z[4]))+1)/int(z[10]),subject_coverage=(abs(int(z[7])-int(z[6]))+1)/int(z[11]),evalue=float(z[8]),bitscore=float(z[9]),exact_full_aa=(mode=='blastp' and qseq[q]==seq[s])))
write(O/'all_alignment_correspondences.tsv',hits)
# Exact identity is checked independently of the similarity-search heuristic.
byhash=defaultdict(list)
for q,m in qm.items():
    if m['representation']=='aa':byhash[hashlib.sha256(qseq[q].encode()).hexdigest()].append(q)
exact=[]
for s,v in seq.items():
    for q in byhash[hashlib.sha256(v.encode()).hexdigest()]:
        if qm[q]['species_taxid']==species(s.split('.')[0]):
            exact.append(dict(**qm[q],irwin_id=s,contamination_score=float(lm[s]['contamination_score'])))
write(O/'exact_sequence_pairs.tsv',exact)
thresholds=[('95id_90cov',95,.9),('90id_80cov',90,.8),('80id_80cov',80,.8),('50id_70cov',50,.7),('30id_50cov',30,.5)]
summary=[];pairs=[]
for cut in [False,True]:
    e=[x for x in exact if not cut or x['contamination_score']>=2]
    summary.append(dict(scope='same_species',score_ge2=cut,threshold='exact_full_aa',identity_pct=100,min_bilateral_coverage=1,literature_proteins=len({x['irwin_id'] for x in e}),catalogue_records=len({x['record_id'] for x in e}),literature_denominator=1286 if cut else 1423))
    for name,pid,cov in thresholds:
        a=[x for x in hits if x['same_species'] and (not cut or x['contamination_score']>=2) and x['identity_pct']>=pid and x['query_coverage']>=cov and x['subject_coverage']>=cov]
        # Union exact comparisons to avoid heuristic omission of a known exact pair.
        li={x['irwin_id'] for x in a}|{x['irwin_id'] for x in e};qs={x['record_id'] for x in a}|{x['record_id'] for x in e}
        summary.append(dict(scope='same_species',score_ge2=cut,threshold=name,identity_pct=pid,min_bilateral_coverage=cov,literature_proteins=len(li),catalogue_records=len(qs),literature_denominator=1286 if cut else 1423))
        if cut:
            for x in a:pairs.append(dict(threshold=name,**x))
write(O/'threshold_summary.tsv',summary);write(O/'threshold_pairs_score_ge2.tsv',pairs)
expected={r['irwin_sequence_id'] for r in read(ROOT/'analysis/revision_pnas/revision_round2_assessment/supplementary_tables_reader_20261006/supporting_data/S3/Irwin_protein_correspondence.tsv')}
observed={r['irwin_id'] for r in exact if r['contamination_score']>=2}
status=dict(author_proteins=len(seq),score_ge2=sum(float(x['contamination_score'])>=2 for x in lit),protein_queries=len(aa),nucleotide_queries=len(nt),exact_proteins=len(observed),exact44_reproduced=observed==expected,missing_expected=sorted(expected-observed),additional_exact=sorted(observed-expected),thresholds_fixed_before_analysis=True,logs_read=False)
(O/'verification.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2));print(json.dumps(summary,indent=2))
