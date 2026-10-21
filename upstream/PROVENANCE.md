# Provenance

Where more than one script exists for a stage, this file says which one ran last.
It is built from the modification times of the scripts as they stand on the HPC,
so it is a record of when each was last edited, not of when each job was
submitted. For the autism arm those two are close together and the ordering is
unambiguous. For the schizophrenia and depression arms the timestamps are all
identical, because the scripts were moved into a final directory in one
operation; there, which script ran is settled by where it was filed instead, and
that is noted below.

The scripts themselves are unmodified. Nothing here has been added to them.

## Autism · PRJNA1022510

The most tangled cohort, and the one where the filenames mislead.

| | Script | Status |
|---|---|---|
| 06-13 21:52 | `download_PRJNA1022510.sh` | ran |
| 06-13 22:46 | `fastqc_PRJNA1022510.sh` | ran |
| 06-15 11:37 | `fix_ids_PRJNA1022510.sh` | ran |
| 06-18 21:51 | `taxonomy_PRJNA1022510.sh` | superseded |
| 06-19 15:31 | `dada2_PRJNA1022510.sh` | superseded |
| 06-21 22:20 | `rebuild2_PRJNA1022510.sh` | ran |
| 06-21 23:04 | `test_orient.sh` | diagnostic |
| 06-21 23:25 | `dada2_final_PRJNA1022510.sh` | **superseded, despite the name** |
| 06-21 23:37 | `test_golayoff.sh` | diagnostic |
| 06-22 00:26 | `dada2_correct_PRJNA1022510.sh` | **last DADA2 run** |
| 06-22 00:46 | `taxonomy_PRJNA1022510_final.sh` | **superseded, despite the name** |
| 06-22 14:49 | `taxonomy_PRJNA1022510_redo.sh` | superseded |
| 06-22 17:48 | `taxonomy_PRJNA1022510_oom_fix.sh` | **last taxonomy run** |

`dada2_final_PRJNA1022510.sh` was replaced an hour after it was written, and
`taxonomy_PRJNA1022510_final.sh` is three versions out of date. Neither name
should be trusted; the timestamps above should.

The interleaving is the interesting part. `rebuild2` and `test_orient` come
before `dada2_final`, and `test_golayoff` comes between `dada2_final` and
`dada2_correct`. The reads were in mixed orientation and it was not known whether
the barcodes needed Golay error correction, so both were tested — `test_orient.sh`
demultiplexes each way, `test_golayoff.sh` counts reads recovered with the
correction on and off — and DADA2 was then rerun against the answer rather than
an assumption. That is why `dada2_correct` exists.

## Autism · PRJNA769228

| | Script | Status |
|---|---|---|
| 06-13 22:17 | `download_PRJNA769228.sh` | ran |
| 06-13 22:46 | `fastqc_PRJNA769228.sh` | ran |
| 06-18 22:22 | `dada2_PRJNA769228.sh` | superseded |
| 06-22 14:46 | `dada2_PRJNA769228_fixed.sh` | superseded |
| 06-26 00:41 | `taxonomy_769228_fixed.sh` | superseded |
| 06-27 00:35 | `cutadapt_test_769228.sh` | diagnostic |
| 06-27 00:49 | `dada2_769228_v4.sh` | **last DADA2 run** |
| 06-27 03:20 | `taxonomy_v4_PRJNA769228.sh` | **last taxonomy run** |

`cutadapt_test_769228.sh` is fourteen minutes before `dada2_769228_v4.sh`. The
primer check is what sent this cohort back to be reprocessed restricted to V4, so
the two `_v4` scripts are the ones behind its feature table and the earlier four
are not.

## Autism · PRJNA453621 and PRJNA529598

One script per stage. Nothing to resolve.

| | | |
|---|---|---|
| 06-13 21:33 / 22:17 | `download_PRJNA453621.sh`, `download_PRJNA529598.sh` | ran |
| 06-13 22:46 | `fastqc_*` | ran |
| 06-18 22:22 / 23:16 | `taxonomy_PRJNA453621.sh`, `taxonomy_PRJNA529598.sh` | ran |
| 06-18 23:12 / 23:23 | `dada2_PRJNA529598.sh`, `dada2_PRJNA453621.sh` | ran |

## Autism · pooled stages and reference

| | Script | Status |
|---|---|---|
| 06-13 23:20 | `multiqc_microbiome.sh` | ran |
| 06-14 04:08 | `reference/download_silva.sh` | ran |
| 06-14 04:15 | `reference/train_silva.sh` | ran — note this is the one stage that used QIIME 2 2025.4 |
| 07-01 15:54 | `qiime_phase1.sh` | ran |
| 07-01 16:27 | `qiime_phase2.sh` | ran |
| 07-06 05:07 | `qiime_sepp.sh` | superseded |
| 07-06 09:31 | `qiime_phase2_sepp.sh` | **last per-cohort SEPP run** |
| 07-14 02:01 | `picrust2_run.sh` | ran |
| 07-15 01:01 | `sepp_pooled.sh` | ran — SEPP on the pooled table, a separate stage |

The top-level README names `qiime_sepp.sh` in its pipeline listing.
`qiime_phase2_sepp.sh` is four hours later and is the one that ran last;
`sepp_pooled.sh` is not a replacement for either, it inserts the pooled table.

## Schizophrenia · PRJNA880407

Every script in this arm carries the same timestamp, `2026-07-29 23:56`, from
being moved into the final directory together. The timestamps cannot order them.
The numbering can, and the directory does: these all sat in `scripts/`, the
working set, while anything replaced had been moved to `ARCHIVE_superseded/`.

| Script | Status |
|---|---|
| `00_setup_880407.sh` | ran — directory and metadata setup, recovered from the project directory it was written in |
| `01_download_880407.sh` | ran |
| `01b_redownload_fullq.sh` | ran — a second retrieval of the same deposit at full quality, which is why it carries the `b` |
| `02_fastqc_880407.sh` | ran |
| `03_dada2_880407.sh` | paired-end attempt |
| `03b_dada2_single.sh` | forward read only |
| `03c_dada2_rev_v4.sh` | **the one that produced the table** |
| `04_taxonomy_880407.sh` | ran |
| `05_diversity_880407.sh` | ran |
| `06_ancombc_880407.sh` | ran |
| `07_sepp_880407.sh` | ran |
| `08_sepp_diversity.sh` | ran |

The three DADA2 scripts are a sequence, not alternatives: paired-end merging
failed, the forward read alone was tried, and it was finally the reverse read —
reverse-complemented, because R2 is what carried V4 in forward orientation — that
was denoised. A fourth script, `ARCHIVE_superseded/ancombc_v1_holm/06_ancombc_880407.sh`,
is byte-identical to the one in `scripts/` and is not included here.

## Depression · PRJNA591924 and PRJNA776170

Same situation: one timestamp for the whole arm. Resolved by filing. Three
scripts were in `Microbiome_MDD/scripts/` and the rest in
`ARCHIVE_superseded/SCRIPTS_ARCHIVE/`.

| Script | Status |
|---|---|
| `download_PRJNA591924.sh`, `download_PRJNA776170.sh` | ran |
| `fastqc_PRJNA591924.sh`, `fastqc_PRJNA776170.sh`, `multiqc_PRJNA591924.sh` | ran |
| `cutadapt_PRJNA776170.sh` | ran |
| `dada2_PRJNA591924.sh`, `dada2_v4_PRJNA776170.sh` | ran |
| `taxonomy_PRJNA591924.sh`, `taxonomy_v4_PRJNA776170.sh` | ran |
| `sepp_fixed.sh` | **kept in `scripts/` — the SEPP run that stands** |
| `sepp_debug.sh`, `sepp_only.sh`, `sepp_phylo_mdd.sh` | superseded, archived at the time |
| `run_picrust2.sh` | ran |

The SEPP attempts are the only ambiguity in this arm, and it was already settled:
`sepp_fixed.sh` is the one that was kept in the working directory, and the other
three had been moved aside before any of this.
