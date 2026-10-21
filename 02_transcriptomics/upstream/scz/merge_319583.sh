#!/bin/bash
#SBATCH --job-name=merge_319583
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err
module purge
module load bear-apps/2022a
module load SciPy-bundle/2022.05-foss-2022a
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA319583 || exit 1
mkdir -p counts logs

python3 << 'PYEOF'
import pandas as pd, glob, os

files = sorted(glob.glob("counts/perbam/*.txt"))
print(f"merging {len(files)} per-bam count files...")

merged = None
for i, f in enumerate(files):
    sample = os.path.basename(f).replace(".txt", "")
    # featureCounts: skip the '# Program' line (header=1); cols = Geneid..Length, then count (last col)
    df = pd.read_csv(f, sep="\t", comment="#", usecols=[0, 6], names=["Geneid", sample], header=0)
    if merged is None:
        merged = df
    else:
        merged = merged.merge(df, on="Geneid", how="outer")
    if (i+1) % 50 == 0:
        print(f"  {i+1}/{len(files)}")

merged = merged.set_index("Geneid")
merged.to_csv("counts/PRJNA319583_gene_counts_matrix.tsv", sep="\t")
print(f"DONE: matrix shape = {merged.shape[0]} genes x {merged.shape[1]} samples")
print(f"any NaN? {merged.isna().any().any()}")
PYEOF
echo "=== merge done: $(date) ==="
