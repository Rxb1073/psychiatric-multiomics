#!/bin/bash
set -euo pipefail

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome
DEST="$BASE/PRJNA880407"
PRJ=PRJNA880407

mkdir -p "$DEST"/raw "$DEST"/metadata "$DEST"/fastqc "$DEST"/qiime2_v4 "$DEST"/logs

echo ">> fetching ENA run table for $PRJ ..."
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PRJ}&result=read_run&fields=run_accession,sample_accession,sample_title,sample_alias,fastq_ftp,fastq_md5,library_strategy,library_layout,read_count&format=tsv" -o "$DEST/metadata/ena_file_info.tsv"

n=$(($(wc -l < "$DEST/metadata/ena_file_info.tsv") - 1))
echo ">> runs found: $n   (expected 100)"

echo ""
echo ">> library strategy / layout:"
tail -n +2 "$DEST/metadata/ena_file_info.tsv" | cut -f7,8 | sort | uniq -c

echo ""
echo ">> group counts from sample_alias:"
tail -n +2 "$DEST/metadata/ena_file_info.tsv" | cut -f4 | sed -E 's/-?[0-9]+//g' | sort | uniq -c

echo ""
echo ">> building QIIME2 metadata"
awk -F'\t' 'FNR==1 { print "sample-id\tdiagnosis\tsubject_code"; next } { run=$1; alias=$4; grp="UNKNOWN"; if (alias ~ /^SZ/) grp="SZ"; else if (alias ~ /^C/) grp="HC"; print run "\t" grp "\t" alias }' "$DEST/metadata/ena_file_info.tsv" > "$DEST/metadata/metadata_qiime2.tsv"

echo ">> metadata group counts:"
tail -n +2 "$DEST/metadata/metadata_qiime2.tsv" | cut -f2 | sort | uniq -c

echo ""
echo ">> UNKNOWN labels (should be none):"
awk -F'\t' 'NR>1 && $2=="UNKNOWN"' "$DEST/metadata/metadata_qiime2.tsv" || true

echo ""
echo "Folder ready: $DEST"
