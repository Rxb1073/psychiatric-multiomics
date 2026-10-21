# Transcriptomics upstream processing

SLURM job scripts that take each RNA-seq deposit from raw FASTQ to a gene count
matrix. Everything downstream — differential expression, meta-analysis,
set-level testing, deconvolution — is in the R scripts one directory up.

These are the scripts as they ran on the University of Birmingham BlueBear HPC,
unedited. Paths, SLURM resource requests and `module load` lines are left in
place: they are the record of what was actually executed, and the resource
requests and job-array ranges are often the only surviving note of how many
samples a stage handled. They will not run elsewhere without the paths being
changed.

## Stage order

```
shared/build_star_index.sh       STAR index, built once and reused
download_<accession>.sh          ENA REST API, fetched with wget, MD5 verified
fastqc_<accession>.sh            per-run FastQC
trimgalore_*.sh / trim_*.sh      adapter and quality trimming
fastqc_posttrim_*.sh             FastQC repeated after trimming
star_<accession>.sh              STAR alignment
featurecounts_*.sh / fc_*.sh     featureCounts, strandedness set per deposit
merge_matrices.sh                per-run counts into one matrix per cohort
```

All deposits were aligned against the same reference: Ensembl GRCh38 release
109, primary assembly, with the matching GTF for counting. One index, built by
`shared/build_star_index.sh`, served every arm, so a difference between deposits
is not a difference in annotation.

`shared/inspect_pheno.sh` pulls each GEO series matrix and prints its phenotype
columns. It is how the donor-level covariates were read before any sample sheet
was written, and it is the step behind the donor counts in the top-level dataset
table.

Strandedness was determined per deposit rather than assumed, and the resulting
`-s` setting for every deposit is recorded in the dataset table in the top-level
README. Getting this wrong silently halves the counts, so it is set explicitly in
each `featurecounts_*.sh` rather than left to a default.

## Directories

| | |
|---|---|
| `shared/` | STAR index, pooled post-trim QC, phenotype staging |
| `asd/` | GSE102741, PRJNA434002, GSE59288 / GSE51264, GSE211154 |
| `scz/` | PRJNA314463, PRJNA319583, PRJNA574470 |
| `mdd/` | GSE101521, GSE102556 |

The schizophrenia deposits appear here under their SRA BioProject accessions
rather than the GEO series accessions used in the top-level README, because that
is how they were downloaded. Two of the three can be matched from evidence inside
the scripts:

| BioProject | GEO series | On what basis |
|---|---|---|
| PRJNA574470 | GSE138082 | Three independent matches. `--array=1-78` and `1-79` against the README's 78 samples; `-s 2`, which is the only one of the three deposits with that strandedness; and `metadata/pe_bams.txt` plus `metadata/se_bams.txt`, matching the mixed PE + SE library layout that only GSE138082 has. |
| PRJNA319583 | GSE80655 | `featurecounts_PRJNA319583.sh` is byte-identical in the schizophrenia and the depression directories. GSE80655 is the only deposit that contributes to both arms, which is the only thing that explains one counting script shared between them. |
| PRJNA314463 | GSE78936 | `-s 0` and, resolving the apparent size mismatch, the series composition: GSE78936 spans three brain regions and includes bipolar donors alongside schizophrenia and control (Hu *et al.* 2016), so a trimming array of `1-82` over the whole deposit is consistent with the 52 samples and 28 BA11 donors the top-level README analyses. The bipolar donors were not carried into any analysis. |

GSE211154 was processed and then excluded from the analysis: it is cerebellar
tissue rather than cortex, and reached only 65.6% median unique alignment. Its
scripts are kept because the exclusion is a result, and the alignment rate behind
it came from this run.

## Reruns and variants

Where more than one script exists for a stage, they are the successive attempts,
kept rather than pruned so that the record matches what ran.

**Which one ran is in `PROVENANCE.md`.** It resolves each stage against the
scripts' modification times on the HPC and says plainly where it cannot — for
PRJNA574470 it cannot, and that is stated rather than guessed.

| Suffix | |
|---|---|
| `_array` | Rewritten as a SLURM job array after the serial version proved too slow |
| `_fixed`, `_perBAM` | Rerun after the previous attempt failed or was incomplete |
| `_150bp` | Reprocessed at the correct read length |
| `*_stragglers`, `retrim_*`, `star_2samples`, `realign_and_fc` | Targeted reruns of the handful of samples that failed in the main job |
| `recover_*` | An attempted rerun that never ran — see `PROVENANCE.md` |
| `namesort_*`, `sort_unsorted_bams` | BAM sort-order repair before featureCounts |

Two things in `PROVENANCE.md` are worth knowing before reading the scripts
themselves. For GSE102741 the name-sorting scripts postdate all four
featureCounts attempts, because sort order was the underlying problem rather than
the counting parameters. And `recover_PRJNA434002_6samples.sh` never ran: it
names six accessions in a series that does not occur in this deposit, left no
output, and PRJNA434002's count matrix holds 41 libraries — the full deposited
total. Nothing is missing from that cohort.

The straggler scripts are the most informative of the set.
`retrim_empty_PRJNA434002.sh` rebuilds its own list of empty outputs by file size
rather than taking a list on trust, and five samples across GSE101521 and
GSE102556 were realigned and counted separately before being merged. Those
samples are in the published counts, and these scripts are how they got there.

## Controls

The schizophrenia arm carries six scripts that are not pipeline stages:
`verify_md5.sh`, `verify_se.sh`, `dup_check.sh`, `sweep_dup.sh`,
`confound_check.sh` and `diag.sh`, with `check_trim.sh` in `shared/`. They verify
the download by checksum, confirm that libraries recorded as single-end really
are, detect and sweep duplicates, and screen for confounds — all before counting.
`PROVENANCE.md` places them in sequence.
