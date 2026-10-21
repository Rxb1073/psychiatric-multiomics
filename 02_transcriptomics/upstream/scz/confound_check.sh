#!/bin/bash
D=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470
SHARED=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/metadata
META=$D/metadata/ena_file_info.tsv
CSV=$(ls $SHARED/*574470*.csv 2>/dev/null | head -1)
module purge; module load bear-apps/2024a; module load Python/3.12.3-GCCcore-13.3.0 2>/dev/null
tail -n +2 "$META" | awk -F'\t' '{n=gsub(/fastq.gz/,"&",$2); print $1"\t"(n==1?"SE":"PE")}' > $D/metadata/layout_by_run.tsv
echo "CSV used: $CSV"
python3 << PYEOF
import csv
from collections import Counter, defaultdict
layout = dict(l.split() for l in open("$D/metadata/layout_by_run.tsv"))
rows = list(csv.DictReader(open("$CSV")))
cols = list(rows[0].keys())
runcol = next((c for c in cols if c.lower() in ('run','run_accession')), None)
dxcol  = next((c for c in cols if 'diagn' in c.lower()), None)
print("run column:", runcol, "| diagnosis column:", dxcol)
print("all columns:", cols)
tab = defaultdict(Counter)
for r in rows:
    tab[layout.get(r.get(runcol,''),'?')][r.get(dxcol,'NA')] += 1
print("\n--- layout x diagnosis ---")
for lay in sorted(tab):
    print(f"{lay}: {dict(tab[lay])}")
PYEOF
