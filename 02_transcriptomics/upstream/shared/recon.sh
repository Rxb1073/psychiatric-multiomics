#!/bin/bash
BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq
for D in PRJNA319583 PRJNA314463 PRJNA574470; do
  echo "================ $D ================"
  DIR=$BASE/$D
  META=$DIR/metadata/ena_file_info.tsv
  if [ -f "$META" ]; then
    RUNS=$(tail -n +2 "$META" | wc -l)
    EXP=$(tail -n +2 "$META" | cut -f2 | tr ';' '\n' | grep -c 'fastq.gz')
    PERRUN=$(tail -n +2 "$META" | head -1 | cut -f2 | tr ';' '\n' | grep -c 'fastq.gz')
    if [ "$PERRUN" -ge 2 ]; then LAYOUT="PAIRED-END"; else LAYOUT="SINGLE-END"; fi
    echo "ENA manifest: $RUNS runs | $EXP expected files | ${PERRUN} file(s)/run -> $LAYOUT"
  else
    echo "!! no ena_file_info.tsv at $META"
    EXP=""
  fi
  GOT=$(find $DIR/fastq_files -name '*.fastq.gz' 2>/dev/null | wc -l)
  TINY=$(find $DIR/fastq_files -name '*.fastq.gz' -size -1M 2>/dev/null | wc -l)
  echo "on disk: $GOT files | <1MB (partial/failed): $TINY | expected: ${EXP:-?}"
  if [ -n "$EXP" ]; then
    if [ "$GOT" -ge "$EXP" ]; then echo ">> COMPLETE"; else echo ">> INCOMPLETE: missing $((EXP-GOT))"; fi
  fi
  R1=$(find $DIR/fastq_files -name '*_1.fastq.gz' 2>/dev/null | wc -l)
  R2=$(find $DIR/fastq_files -name '*_2.fastq.gz' 2>/dev/null | wc -l)
  echo "naming: _1=$R1  _2=$R2"
  LOG=$(ls -t $DIR/logs/dl_*.out 2>/dev/null | head -1)
  if [ -n "$LOG" ]; then
    echo "md5 verdict (from $(basename $LOG)):"
    grep -E 'VERIFIED|corrupted|FAILED' "$LOG" | tail -2
  fi
  F=$(find $DIR/fastq_files -name '*.fastq.gz' 2>/dev/null | head -1)
  if [ -n "$F" ]; then
    echo "length histogram (first 1000 reads of $(basename $F)):"
    zcat "$F" 2>/dev/null | head -4000 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c
  fi
  echo ""
done
