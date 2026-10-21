#!/bin/bash
D=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470
META=$D/metadata/ena_file_info.tsv

echo "--- files per run (2=PE, 1=SE) ---"
tail -n +2 "$META" | awk -F'\t' '{n=gsub(/fastq.gz/,"&",$2); print n}' | sort | uniq -c

echo ""
echo "--- SE runs (single file) ---"
tail -n +2 "$META" | awk -F'\t' '{n=gsub(/fastq.gz/,"&",$2); if(n==1) print $1}'

echo ""
echo "--- missing file(s): in manifest, not on disk ---"
tail -n +2 "$META" | cut -f2 | tr ';' '\n' | grep fastq.gz | sed 's|.*/||' | sort -u > /tmp/exp574.txt
ls $D/fastq_files/*.fastq.gz 2>/dev/null | sed 's|.*/||' | sort -u > /tmp/got574.txt
comm -23 /tmp/exp574.txt /tmp/got574.txt
