#!/bin/bash
D=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470
META=$D/metadata/ena_file_info.tsv
FQ=$D/fastq_files

echo "=== (1) confirm the 10 manifest-SE runs are bare-file-only (not broken PE) ==="
SE_RUNS=$(tail -n +2 "$META" | awk -F'\t' '{n=gsub(/fastq.gz/,"&",$2); if(n==1) print $1}')
for S in $SE_RUNS; do
  bare=$(ls $FQ/${S}.fastq.gz 2>/dev/null | wc -l)
  p1=$(ls $FQ/${S}_1.fastq.gz 2>/dev/null | wc -l)
  p2=$(ls $FQ/${S}_2.fastq.gz 2>/dev/null | wc -l)
  if [ "$bare" -eq 1 ] && [ "$p1" -eq 0 ] && [ "$p2" -eq 0 ]; then
    echo "OK   $S : bare SE"
  else
    echo "FLAG $S : bare=$bare _1=$p1 _2=$p2  <-- not clean SE, investigate"
  fi
done

echo ""
echo "=== (2) read length: SE vs PE branch ==="
SE1=$(echo "$SE_RUNS" | head -1)
PE1=$(ls $FQ/*_1.fastq.gz 2>/dev/null | head -1 | sed 's|.*/||; s|_1.fastq.gz||')
echo "SE  $SE1 : $(zcat $FQ/${SE1}.fastq.gz   | head -400 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | tr '\n' ' ')"
echo "PE  $PE1 (_1): $(zcat $FQ/${PE1}_1.fastq.gz | head -400 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | tr '\n' ' ')"
