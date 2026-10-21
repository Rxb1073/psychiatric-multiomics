#!/bin/bash
#SBATCH --job-name=picrust
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
set -euo pipefail
DATASET=$1
BASE=/rds/projects/e/elhamsak-dissertation-autism/Microbiome
cd "$BASE/$DATASET" || exit 1
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$BASE/conda_envs/picrust2"
export TMPDIR="$BASE/tmp/picrust_${SLURM_JOB_ID}"; mkdir -p "$TMPDIR"

rm -rf picrust2_out        # picrust2 refuses to write into an existing dir
echo "=== $DATASET PICRUSt2 | $(date) ==="

picrust2_pipeline.py \
  -s picrust2_input/seqs/dna-sequences.fasta \
  -i picrust2_input/tab/feature-table.biom \
  -o picrust2_out \
  -p 8 \
  --stratified \
  --verbose

echo "=== NSTI check ==="
zcat picrust2_out/*_marker_predicted_and_nsti.tsv.gz 2>/dev/null | head -1
# mean weighted NSTI (prediction quality; <0.15 good)
python3 - <<'PY'
import gzip, glob
f=glob.glob("picrust2_out/marker_predicted_and_nsti.tsv.gz")+glob.glob("picrust2_out/*nsti*.tsv.gz")
if f:
    import statistics
    vals=[]
    with gzip.open(f[0],'rt') as fh:
        h=fh.readline().split('\t'); 
        idx=[i for i,c in enumerate(h) if 'nsti' in c.lower()]
        if idx:
            for line in fh:
                try: vals.append(float(line.split('\t')[idx[0]]))
                except: pass
    if vals: print(f"mean NSTI = {statistics.mean(vals):.3f}  (n={len(vals)}, >0.15: {sum(v>0.15 for v in vals)})")
PY
echo "=== $DATASET done | $(date) ==="
ls -lh picrust2_out/pathways_out/path_abun_unstrat.tsv.gz
