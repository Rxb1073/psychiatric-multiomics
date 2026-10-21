#!/bin/bash
D=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470/fastq_files
DUP=0; GEN=0
for R1 in $D/*_1.fastq.gz; do
  S=$(basename $R1 _1.fastq.gz)
  R2=$D/${S}_2.fastq.gz
  [ -f "$R2" ] || continue
  A=$(zcat $R1 | head -40000 | awk 'NR%4==2' | md5sum | cut -d' ' -f1)
  B=$(zcat $R2 | head -40000 | awk 'NR%4==2' | md5sum | cut -d' ' -f1)
  if [ "$A" = "$B" ]; then echo "DUPLICATE _2: $S"; DUP=$((DUP+1)); else GEN=$((GEN+1)); fi
done
echo ""
echo ">> genuine PE: $GEN  |  _2-duplicates-_1: $DUP"
