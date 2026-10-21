#!/bin/bash
#SBATCH --job-name=merge_574470
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err
module purge
module load bear-apps/2022a
module load SciPy-bundle/2022.05-foss-2022a
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p counts logs
python3 << 'PYEOF'
import pandas as pd, glob, os
files = sorted(glob.glob("counts/perbam/*.txt"))
print(f"merging {len(files)} per-bam files (PE + SE together)...")
merged = None
for f in files:
    sample = os.path.basename(f).replace(".txt", "")
    df = pd.read_csv(f, sep="\t", comment="#", usecols=[0, 6], names=["Geneid", sample], header=0)
    merged = df if merged is None else merged.merge(df, on="Geneid", how="outer")
merged = merged.set_index("Geneid")
merged.to_csv("counts/PRJNA574470_gene_counts_matrix.tsv", sep="\t")
print(f"DONE: {merged.shape[0]} genes x {merged.shape[1]} samples")
print(f"any NaN? {merged.isna().any().any()}")
PYEOF
echo "=== merge done: $(date) ==="
