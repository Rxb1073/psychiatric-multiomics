#!/bin/bash
D=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470/fastq_files
S=$(ls $D/*_2.fastq.gz 2>/dev/null | head -1 | sed 's|.*/||; s|_2.fastq.gz||')
echo "Testing a run that has BOTH files: $S"
echo ""
echo "=== read IDs + sequences, first 2 reads each ==="
echo "--- ${S}_1 ---"; zcat $D/${S}_1.fastq.gz | head -8
echo "--- ${S}_2 ---"; zcat $D/${S}_2.fastq.gz | head -8
echo ""
echo "=== sequence-only fingerprint (first 10k reads) ==="
A=$(zcat $D/${S}_1.fastq.gz | head -40000 | awk 'NR%4==2' | md5sum | cut -d' ' -f1)
B=$(zcat $D/${S}_2.fastq.gz | head -40000 | awk 'NR%4==2' | md5sum | cut -d' ' -f1)
echo "_1 seq md5: $A"
echo "_2 seq md5: $B"
if [ "$A" = "$B" ]; then
  echo ">> IDENTICAL — _2 duplicates _1. Dataset is SINGLE-END; drop all _2."
else
  echo ">> DIFFERENT — check the read IDs above: same ID base with /1 vs /2 (or '1:' vs '2:') = genuine paired-end."
fi
