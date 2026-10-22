#!/usr/bin/env bash
set -euo pipefail
source ./config.sh
mkdir -p "$OUTDIR/logs"

find "$OUTDIR/jobs" -mindepth 2 -maxdepth 2 -name input.aa.fa > "$OUTDIR/tree_input_fastas.list"

cat > "$OUTDIR/run_one_tree.sh" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
fa="$1"
threads="$2"
jobdir=$(dirname "$fa")
cd "$jobdir"
seq_n=$(grep -c '^>' input.aa.fa || true)
if [[ "$seq_n" -lt 4 ]]; then
  echo "SKIP $jobdir: only $seq_n sequences" > tree.status
  exit 0
fi
mafft --auto --thread "$threads" input.aa.fa > aln.fa 2> mafft.log
trimal -in aln.fa -out aln.trim.fa -automated1 > trimal.log 2>&1
trim_n=$(grep -c '^>' aln.trim.fa || true)
if [[ "$trim_n" -lt 4 ]]; then
  echo "SKIP $jobdir: only $trim_n sequences after trimming" > tree.status
  exit 0
fi
iqtree2 -s aln.trim.fa -m MFP -B 1000 --alrt 1000 -T "$threads" -pre iqtree > iqtree.stdout 2> iqtree.stderr
echo "DONE $jobdir" > tree.status
EOS
chmod +x "$OUTDIR/run_one_tree.sh"

if command -v parallel >/dev/null 2>&1; then
  parallel -j "$PARALLEL_JOBS" "$OUTDIR/run_one_tree.sh" {} "$THREADS" :::: "$OUTDIR/tree_input_fastas.list"
else
  while read -r fa; do
    "$OUTDIR/run_one_tree.sh" "$fa" "$THREADS"
  done < "$OUTDIR/tree_input_fastas.list"
fi
