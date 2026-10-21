#!/bin/bash
#SBATCH --job-name=merge_mdd
#SBATCH --account=elhamsak-2-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/merge_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/merge_%j.err
module purge; module load bear-apps/2022a; module load SciPy-bundle/2022.05-foss-2022a
python3 << 'PY'
import glob, os, pandas as pd
A2='/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics'
for D in ['GSE102556','GSE101521']:
    perbam=f'{A2}/{D}/counts/perbam'
    files=sorted(glob.glob(f'{perbam}/*.txt'))
    merged=None
    for f in files:
        sample=os.path.basename(f).replace('.txt','')
        df=pd.read_csv(f, sep='\t', comment='#')
        sub=df[['Geneid', df.columns[-1]]].rename(columns={df.columns[-1]: sample})
        merged=sub if merged is None else merged.merge(sub, on='Geneid', how='outer')
    os.makedirs(f'{A2}/{D}/counts', exist_ok=True)
    out=f'{A2}/{D}/counts/{D}_gene_counts_matrix.tsv'
    merged.to_csv(out, sep='\t', index=False)
    print(f'{D}: {merged.shape[0]} genes x {merged.shape[1]-1} samples -> {out}')
PY
