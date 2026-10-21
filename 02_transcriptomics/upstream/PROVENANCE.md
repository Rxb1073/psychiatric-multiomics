# Provenance

Where more than one script exists for a stage, this file says which one ran last.
It is built from the modification times of the scripts as they stand on the HPC,
so it is a record of when each was last edited, not of when each job was
submitted. The autism and depression arms order cleanly. The schizophrenia arm's
scripts all carry one timestamp, from being copied into a final directory
together; that is noted below.

The scripts themselves are unmodified. Nothing here has been added to them.

## Shared

| | Script | |
|---|---|---|
| 06-15 08:05 | `shared/multiqc_RNAseq.sh` | ran |
| 06-15 08:07 | `shared/build_star_index.sh` | ran — Ensembl GRCh38 release 109, built once, used by every arm |
| 06-19 22:27 | `shared/fastqc_posttrim.sh` | ran |
| 06-28 21:01 | `shared/inspect_pheno.sh` | ran — GEO series matrices read for donor-level covariates |

## Autism · GSE102741

| | Script | Status |
|---|---|---|
| 06-14 04:19 | `download_GSE102741.sh` | ran |
| 06-14 16:19 | `fastqc_GSE102741.sh` | ran |
| 06-18 23:36 | `trimgalore_GSE102741.sh` | ran |
| 06-19 21:04 | `star_GSE102741.sh` | ran |
| 06-22 12:09 | `featurecounts_GSE102741.sh` | superseded |
| 06-26 01:12 | `featurecounts_GSE102741_fixed.sh` | superseded |
| 06-26 21:33 | `featurecounts_GSE102741_array.sh` | superseded |
| 06-27 13:46 | `featurecounts_GSE102741_perBAM.sh` | **last featureCounts script written** |
| 06-27 17:02 | `namesort_test.sh` | diagnostic |
| 06-27 17:54 | `namesort_remaining.sh` | ran, and postdates every featureCounts script above |

Read this one carefully. Four featureCounts attempts, and then the two
name-sorting scripts come *after* all of them. featureCounts needs name-sorted
BAMs for paired-end counting; the sort order was the underlying problem, and
`namesort_remaining.sh` is the fix. `featurecounts_GSE102741_perBAM.sh` is the
counting script that stands, but it had to be run again once the BAMs were
sorted, and that rerun left no new script behind. The sequence is
`namesort_remaining.sh` and then `featurecounts_GSE102741_perBAM.sh`, not the
order the timestamps are printed in.

## Autism · PRJNA434002

| | Script | Status |
|---|---|---|
| 06-14 04:39 | `download_PRJNA434002.sh` | ran |
| 06-14 16:19 | `fastqc_PRJNA434002.sh` | ran |
| 06-20 13:33 | `trimgalore_PRJNA434002.sh` | superseded |
| 06-20 17:39 | `trimgalore_PRJNA434002_array.sh` | **last trimming script** |
| 06-21 10:51 | `star_PRJNA434002.sh` | superseded |
| 06-22 01:03 | `retrim_empty_PRJNA434002.sh` | ran — rebuilds its own list of empty outputs by file size |
| 06-22 14:20 | `star_PRJNA434002_150bp.sh` | **main alignment** |
| 06-26 00:32 | `star_2samples.sh` | ran — two samples that failed the main job |
| 06-26 01:20 | `featurecounts_PRJNA434002_150bp.sh` | **main counting** |
| 06-26 11:03 | `realign_and_fc.sh` | ran — realign plus count, same day |
| 07-21 08:20 | `recover_PRJNA434002_6samples.sh` | **never ran — see below** |

`recover_PRJNA434002_6samples.sh` is dated 21 July, three weeks after the rest of
this cohort, but it never ran. It names six accessions in the `SRR9262` series,
while every library in PRJNA434002 is `SRR9292` — including those named in
`star_2samples.sh` and `realign_and_fc.sh` — and `SRR9262` appears nowhere else
in the project. It also points at a different project directory from the rest of
this cohort's scripts. It left no output, and the cohort's count matrix
holds 41 libraries, matching the deposited total. Nothing is missing from this
cohort; the script is a record of an attempt that could not have succeeded.

The cohort took three completed rounds: the array rewrite of trimming, the
150 bp realignment, and the two-sample and empty-output repairs, each of which
did write results.

## Autism · GSE59288 / GSE51264 and GSE211154

One script per stage, written in the same sittings. Nothing to resolve.

| | | |
|---|---|---|
| 06-14 04:19 | `download_GSE59288_GSE51264.sh`, `download_GSE211154.sh` | ran |
| 06-14 16:54 | `fastqc_*` | ran |
| 06-15 09:02 | `trimgalore_GSE59288.sh`, `trimgalore_GSE211154.sh` | ran |
| 06-19 21:04 | `star_*` | ran |
| 06-21 11:03 | `featurecounts_*` | ran |

GSE211154 was processed and then excluded: cerebellar tissue rather than cortex,
and 65.6% median unique alignment. Its scripts are kept because the exclusion is
a result and the alignment rate behind it came from this run.

## Depression · GSE101521 and GSE102556

A single clean pass with stragglers at the end. No superseded scripts in this arm.

| | Script | |
|---|---|---|
| 07-10 23:02 | `download_GSE101521.sh`, `download_GSE102556.sh` | ran |
| 07-11 12:32 | `fastqc_GSE101521.sh`, `fastqc_GSE102556.sh` | ran |
| 07-11 12:54 | `multiqc_MDD_rna.sh` | ran |
| 07-11 13:09 | `trim_GSE101521.sh`, `trim_GSE102556.sh` | ran |
| 07-11 21:49 | `star_GSE101521.sh`, `star_GSE102556.sh` | ran |
| 07-11 23:15 | `fc_GSE101521.sh`, `fc_GSE102556.sh` | ran |
| 07-12 00:10 | `sort_unsorted_bams.sh` | ran |
| 07-12 00:39 | `realign_5_stragglers.sh` | ran — five samples across both deposits |
| 07-12 01:07 | `fc_stragglers.sh` | ran |
| 07-12 01:55 | `merge_matrices.sh` | ran — final count matrix per cohort |

Fourteen hours end to end, and the straggler handling is explicit rather than
silent: five samples were realigned and counted separately, then merged with the
rest. Those five are in the published counts.

## Schizophrenia · PRJNA314463, PRJNA319583, PRJNA574470

Recovered from the project directory these were originally written in, so these
timestamps are real edit times rather than the single copy-date that the
`DISSERTATION_final` copies carry.

| | Script | Status |
|---|---|---|
| 06-26 15:10 | `download_rnaseq.sh` | ran — shared retrieval driver |
| 06-26 22:47 | `download_PRJNA314463.sh`, `download_PRJNA319583.sh`, `download_PRJNA574470.sh` | ran |
| 06-27 17:06 | `recon.sh` | diagnostic — deposit reconnaissance |
| 06-27 17:10 | `fastqc_PRJNA314463.sh`, `fastqc_PRJNA319583.sh` | ran |
| 06-27 17:12 | `diag.sh` | diagnostic |
| 06-27 17:17 | `dup_check.sh` | **control — duplicate detection** |
| 06-27 17:19 | `sweep_dup.sh` | ran — acts on the duplicate check |
| 06-27 17:21 | `verify_se.sh` | **control — confirms single-end libraries are single-end** |
| 06-27 17:37 | `confound_check.sh` | **control — confound screen** |
| 06-27 17:52 | `verify_md5.sh` | **control — checksum verification of the download** |
| 06-27 18:36 | `fastqc_PRJNA574470.sh` | ran |
| 06-27 21:08 | `trim_array_PRJNA314463.sh`, `trim_array_PRJNA319583.sh`, `trim_array_PRJNA574470.sh` | ran |
| 06-28 00:37 | `check_trim.sh` | control — trimming verification |
| 06-28 01:04 | `fastqc_posttrim_*` | ran |
| 06-28 02:57 | `multiqc_posttrim.sh` | ran — in `shared/` |
| 06-28 03:59 | `star_PRJNA314463.sh`, `star_PRJNA319583.sh`, `star_PRJNA574470.sh` | ran |
| 06-28 06:08 | `featurecounts_PRJNA314463.sh`, `featurecounts_PRJNA574470.sh` | PRJNA314463 ran; PRJNA574470 superseded |
| 06-28 15:26 | `featurecounts_PRJNA319583.sh` | ran |
| 06-28 17:49 | `merge_319583.sh` | ran |
| 06-29 01:07 | `featurecounts_array_PRJNA574470.sh` | superseded by the line below |
| 06-29 16:02 | `fc_heavy5.sh` | **last counting step for PRJNA574470** |
| 06-29 17:30 | `merge_574470.sh` | ran |

PRJNA574470's counting is a sequence, not a choice: the serial script, then the
job array, then `fc_heavy5.sh` for the five heaviest libraries, then the merge.
Nothing here is an alternative to anything else.

Six of these scripts are controls rather than pipeline stages, and they are the
reason this arm's deposit handling can be checked rather than taken on trust:
the download was verified by checksum, the single-end libraries were confirmed to
be single-end before being treated as such, duplicates were detected and swept,
and a confound screen was run before any counting. They are kept for that reason.
