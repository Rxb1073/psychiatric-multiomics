# Microbiome upstream processing

SLURM job scripts that take each 16S deposit from raw FASTQ to a feature table,
taxonomy, and an inserted phylogeny. Everything downstream of this — batch
correction, diversity, differential abundance, the function bridge — is in the
R and Python scripts one directory up.

These are the scripts as they ran on the University of Birmingham BlueBear HPC,
unedited. Paths, SLURM resource requests and `module load` lines are left in
place: they are the record of what was actually executed, and the resource
requests are often the only surviving note of why a stage was split or rerun.
They will not run elsewhere without the paths being changed.

## Stage order

```
reference/download_silva.sh      SILVA 138 reference, and SILVA 128 for SEPP
reference/train_silva.sh         V4-trimmed naive Bayes classifier

download_<accession>.sh          ENA REST API, fetched with wget, MD5 verified
fastqc_<accession>.sh            per-run FastQC
multiqc_*.sh                     aggregate
cutadapt_*.sh                    primer removal, where primers were still present
dada2_<accession>.sh             denoise, merge, chimera removal, per cohort
taxonomy_<accession>.sh          classify against the trained SILVA classifier
qiime_phase1.sh                  de novo tree, table summary, alpha rarefaction
qiime_phase2.sh                  genus collapse, prevalence filter
qiime_sepp.sh / sepp_*.sh        fragment insertion into the SILVA 128 tree
picrust2_run.sh                  predicted gene content, MetaCyc and KEGG
```

DADA2 runs per cohort, not pooled, because read lengths and primer sets differ
between deposits; the cohorts are only brought together after taxonomy, at genus
level.

## Directories

| | |
|---|---|
| `reference/` | Shared reference data and classifier training |
| `asd/` | PRJNA453621, PRJNA529598, PRJNA769228, PRJNA1022510, and the pooled autism stages |
| `scz/` | PRJNA880407, numbered `01b`–`08` |
| `mdd/` | PRJNA591924, PRJNA776170 |

The three arms were processed in separate project directories at different
points in the project, which is why the naming conventions differ between them:
the autism arm is `<stage>_<accession>.sh`, the schizophrenia arm is numbered by
stage, and the depression arm sits between the two. The pipeline steps are the
same in each.

## Reruns and variants

Where more than one script exists for a stage, they are the successive attempts,
kept rather than pruned so that the record matches what ran.

**Which one ran is in `PROVENANCE.md`, and the filenames are not a reliable
guide.** For PRJNA1022510, `dada2_final_PRJNA1022510.sh` was replaced an hour
after it was written and `taxonomy_PRJNA1022510_final.sh` is three versions out
of date. `PROVENANCE.md` resolves every stage of every cohort against the
scripts' modification times on the HPC, and says plainly where it cannot.

The suffixes still carry their intended meaning, as a description of why a rerun
happened rather than as a claim about ordering:

| Suffix | |
|---|---|
| `_fixed`, `_correct`, `_final`, `_redo` | Rerun after the previous attempt was wrong or incomplete |
| `_oom_fix` | Rerun after an out-of-memory kill, with the resource request raised |
| `_v4` | Reprocessed restricted to the V4 region |
| `_array` | Rewritten as a SLURM job array |
| `test_*`, `*_test_*` | Diagnostic run, not a pipeline stage |

The diagnostics are worth reading rather than skipping. `test_orient.sh` and
`test_golayoff.sh` exist because PRJNA1022510's reads were in mixed orientation
and it was not known whether the barcodes needed Golay error correction;
`test_orient.sh` demultiplexes each way and `test_golayoff.sh` counts reads
recovered with the correction on and off. DADA2 was then rerun against the
answer. Likewise `cutadapt_test_769228.sh` runs fourteen minutes before
`dada2_769228_v4.sh`: the primer check is what sent that cohort back to be
reprocessed on V4.

For PRJNA880407 the three DADA2 scripts are a sequence, not alternatives.
`03_dada2_880407.sh` is paired-end; `03b_dada2_single.sh` falls back to the
forward read alone; `03c_dada2_rev_v4.sh` reverse-complements every reverse read
and denoises those, because it was R2 that carried V4 in forward orientation.
The third is the one behind the table.

## Schizophrenia setup

`00_setup_880407.sh` and `01_download_880407.sh` precede the numbered stages and
were recovered from the project directory they were written in.
`01b_redownload_fullq.sh` is a second retrieval of the same deposit at full
quality, which is why it carries the `b`.
