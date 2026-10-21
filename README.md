# Multi-omic reanalysis of autism, schizophrenia and major depression

Analysis code for an MSc Bioinformatics dissertation. Twenty-three open-access
deposits were retrieved and twenty-one analysed across three molecular layers —
gut microbial 16S rRNA amplicon sequencing, post-mortem cortical transcriptomics,
and post-mortem cortical DNA methylation and histone modification — with each
disorder processed through the same pipeline so that differences between them
reflect biology rather than the analytical choices of separate published studies.

**No primary data are included in this repository.** All deposits are public and
listed below with their accessions. Sample sheets are excluded because they carry
donor-level covariates. One summary results table is included —
`01_microbiome/differential_abundance_both_methods.csv`, genus-level differential
abundance with no donor-level information — because the concordance requirement
between ANCOM-BC2 and MaAsLin2 cannot be checked without it.

---

## What was found

Three disorder-specific brain signatures, none extending beyond its own disorder:
GABAergic interneuron markers reduced in schizophrenia, immediate-early genes in
depression, and CHD8-bound promoters displaced in autism. One convergence, in the
gut, between schizophrenia and depression.

That convergence is three genera — *Faecalibacterium*, *Fusicatenibacter* and
*Lachnospiraceae* ND3007 group — recurring between the schizophrenia and
depression concordant sets, 3 of 4 and 26 against a shared tested pool of 116
(hypergeometric P = 0.035), all three in the same direction. The pool is built by
intersecting one arm's MaAsLin2-tested genera with the other's ANCOM-BC2-tested
genera, so it is sensitive to which method's tested set is taken per arm; the
pool size is recorded alongside the p-value in
`12_cross_arm_microbiome/cross_arm_overlap.csv` for that reason. The comparable
autism pairings give no overlap worth testing: 2 of a pool of 119 against
schizophrenia (P = 0.28) and none at all against depression.

Three analytical controls each converted a nominally significant result into a
null: resolving deposits to independent donors, estimating rather than assuming
inter-gene correlation, and comparing gene sets against abundance-matched rather
than transcriptome-wide comparators.

---

## Datasets

| Disorder | Accession | Tissue / region | Assay | Deposited | Analysed | Case / control | Role | Notes |
|---|---|---|---|---|---|---|---|---|
| **GUT MICROBIOME** — faecal 16S rRNA amplicon sequencing · 7 BioProjects · 936 samples downloaded, 919 analysed ||||||||| 
| Autism | PRJNA453621 | Faecal | 16S V4 | 286 | 286 | 143 / 143 | primary (pooled) | 2,490 ASVs; 349 genera |
| | PRJNA529598 | Faecal | 16S V4 | 38 | 38 | 18 / 20 | primary (pooled) | 908 ASVs; 204 genera |
| | PRJNA769228 | Faecal | 16S V4 | 198 | 198 | 138 / 60 | primary (pooled) | 3,466 ASVs; 254 genera |
| | PRJNA1022510 | Faecal | 16S V4 | 139 | 139 | 97 / 42 | primary (pooled) | 2,740 ASVs; 290 genera |
| | **Pooled** | Faecal | 16S V4 | 661 | 661 | 396 / 265 | primary | 9,604 ASVs; 429 genera (166 shared) |
| Schizophrenia | PRJNA880407 | Faecal | 16S V4 | 100 | 100 | 50 / 50 | primary | 5,885 ASVs; 540 genera |
| Depression | PRJNA591924 | Faecal | 16S V4 | 90 | 90 | 43 / 47 | primary (pooled) | 1,950 ASVs; 254 genera |
| | PRJNA776170 | Faecal | 16S V4 | 85 | 68 | 24 / 44 | primary (pooled) | 2,473 ASVs; 345 genera; 17 samples excluded |
| | **Pooled** | Faecal | 16S V4 | 175 | 158 | 67 / 91 | primary | 3,765 ASVs; 400 genera (199 shared) |
| **TRANSCRIPTOMICS** — post-mortem brain bulk RNA-seq · 9 deposits · GSE80655 processed once, subset twice ||||||||| 
| Autism | GSE102741 | DLPFC | PE | 52 samples | 52 donors | 13 / 39 | primary | −s 2 |
| | PRJNA434002 | ACC + PFC | PE | 41 samples | 31 donors | 15 / 16 | primary | −s 2; 10 donors sampled at both regions (within-donor ρ = 0.508) |
| | GSE59288 / GSE51264 | PFC, SFG | SE | 151 runs | 72 humans | 34 / 38 | comparison | −s 0; 40 *M. mulatta* and 39 *P. troglodytes* removed; cases and controls from different series |
| | GSE211154 | Cerebellum | SE | 39 samples | 39 | 20 / 19 | **excluded** | −s 1; laser-captured Purkinje neurons; 65.6% median unique alignment |
| Schizophrenia | GSE80655 | DLPFC | poly(A), PE | 57 runs | 48 donors | 24 / 24 | primary (frontal) | −s 0; 9 libraries across two runs each |
| | GSE78936 | BA11 | poly(A), PE | 52 samples | 40 donors (BA11 28) | 16 / 12 | primary (frontal) | −s 0; 12 donors at BA24 and BA9; BA11 subset never duplicated; deposited as PRJNA314463; the series spans three regions and includes bipolar donors (Hu *et al.* 2016), who are not analysed here, the likeliest explanation for the trimming array covering 82 runs against 52 samples; the 30-run difference has not been itemised |
| | GSE138082 | Hippocampus | PE + SE | 78 samples | 26 donors | 13 / 13 | separate region | −s 2; every donor sampled at CA1, CA3 and DG |
| Depression | GSE101521 | BA9 | rRNA-depleted, PE | 86 runs | 59 donors | 30 / 29 | primary (pooled) | −s 2; two platforms, no donor overlap; 27 microRNA libraries excluded |
| | GSE102556 | dlPFC BA8/9 | rRNA-depleted, PE | 343 runs | 48 donors | 26 / 22 | primary (pooled) | −s 1; runs collapse to 263 samples, each donor sampled at 4–6 regions; DLPFC subset gives one sample per donor |
| | GSE80655 | DLPFC | poly(A), PE | 61 runs | 47 donors | 23 / 24 | sensitivity | −s 0; 14 libraries across two runs each; control donors shared with the schizophrenia arm |
| **DNA METHYLATION AND HISTONE MODIFICATION** — post-mortem brain · 7 datasets ||||||||| 
| Autism | GSE61492 | Neural progenitors | CHD8 ChIP-seq | 5 | 5 | — | binding map | 3 CHD8 preparations, 1 CHD7 specificity control, 1 input; no diagnosis contrast |
| | GSE57369 | Mixed / NSC | CHD8 ChIP-seq | — | — | — | **excluded** | No input control in the series; no-input peak calling reached 7.9% and 12.7% precision against the authors' published calls, below a 70% criterion fixed in advance |
| | GSE131706 | Subventricular zone | 450K array | 34 | 34 | 17 / 17 | primary | Analysed separately from GSE80017 (different tissue) |
| | GSE80017 | Prefrontal cortex | 450K array | 18 | 16 | 8 / 8 | primary | Two arrays removed for incomplete two-colour measurement |
| | **Combined** | SVZ + PFC | 450K array | — | 50 | 25 / 25 | reanalysis | Pooled reanalysis of GSE131706 and GSE80017 |
| Schizophrenia | GSE74193 | Prefrontal cortex | 450K array | 675 arrays | 480 subjects | 185 / 295 | primary | 675 arrays resolve to 515 subjects; prenatal cluster then removed |
| | GSE215991 | Prefrontal cortex | H3K9me2 ChIP-seq | 30 | 30 | 15 / 15 | primary | Tested genome-wide and promoter-restricted |
| Depression | GSE41826 | Frontal cortex, NeuN+ | 450K array | 145 | 58 | 29 / 29 | primary | Sorted neuronal fraction only; series also contains glial, bulk and mixed-fraction arrays |

Two deposits were excluded before analysis: **GSE211154** (cerebellar tissue, not
cortex) and **GSE57369** (no input control deposited; peak calling reached 7.9%
and 12.7% precision against the authors' published calls, below a 70% criterion
fixed in advance). **GSE59288** and **GSE51264** are two accessions forming one
composite dataset. **GSE80655** contributes to the schizophrenia frontal analysis
and serves as a depression sensitivity cohort.

---

## Pipeline order

### Microbiome

Raw FASTQ to feature table, per cohort, in `01_microbiome/upstream/` — download,
FastQC, primer removal, DADA2, taxonomy, phylogeny. That directory has its own
README for stage order and the three arms' differing naming conventions, and a
`PROVENANCE.md` that names which script ran where a stage was attempted more than
once. The filenames are not a reliable guide to that and `PROVENANCE.md` is.

```
download_<accession>.sh  ENA REST API, wget, MD5 verified
dada2_<accession>.sh     denoise, merge, chimera removal, per cohort
taxonomy_<accession>.sh  classify against a V4-trimmed SILVA 138 classifier
qiime_phase1.sh          de novo tree, table summary, alpha rarefaction
qiime_phase2.sh          genus collapse, prevalence filter
qiime_phase2_sepp.sh     fragment insertion into the SILVA 128 reference tree
stage3_mmuphin.R         batch correction, diagnosis retained as covariate
stage4_diversity.R       alpha and beta diversity, PERMANOVA with PERMDISP
stage5_diffabund.R       ANCOM-BC2 and MaAsLin2; concordance required
lodo_check.R             leave-one-dataset-out
picrust2_run.sh          predicted gene content, MetaCyc and KEGG
bridge_*.py              taxonomy-to-function bridge, one per arm
```

### Transcriptomics

Raw FASTQ to count matrix, per deposit, in `02_transcriptomics/upstream/`, with
its own README and `PROVENANCE.md`. All deposits were aligned against one STAR
index built from Ensembl GRCh38 release 109, so a difference between deposits is
not a difference in annotation.

```
build_star_index.sh      GRCh38 release 109, built once and reused
trimgalore_*.sh          adapter and quality trimming
trim_array_*.sh          the same as a SLURM job array, where the cohort was large
star_*.sh                STAR alignment
featurecounts_*.sh       featureCounts, strandedness set explicitly per deposit
fc_*.sh                  the same, depression arm naming
merge_matrices.sh        per-run counts into one matrix per cohort
asd_child_confoundfree.R autism primary, two confound-free cohorts
asd_metaanalysis.R       fixed-effect inverse-variance meta-analysis
asd_setlevel_camera.R    marker-set testing, correlation estimated
scz_pooled_frontal.R     schizophrenia frontal pool, 76 donors
scz_camera_composition_adjusted.R   marker sets with composition covariates
scz_composition_donorcorrect.R      deconvolution, two reference signatures
mdd_primary_2cohort.R    depression primary, two cohorts
mdd_filtered_complete.R  module scores, robustness, manner-of-death contrasts
mdd_ieg_final_checks.R   within-cohort permutation null
```

### Epigenomics
```
02_build_references.sh      bowtie2 indices and blacklists, hg19 and hg38
job1_GSE57369_chipseq.sh    CHD8 ChIP-seq, the series later excluded
job2_GSE61492_chipseq.sh    CHD8 peak calling (MACS2), three-antibody consensus,
                            CHD7 specificity control
job2b_GSE61492_bigwigs.sh   coverage tracks
job3_GSE215991_chipseq.sh   H3K9me2 alignment and peak calling
job_csaw_GSE215991.sh       csaw sliding windows, genome-wide (see Known gaps)
GSE215991_promoter_test.R   the same windows restricted to promoters
GSE215991_sex_control.R     sex as a negative control
GSE215991_rizavi_genes.R    the source publication's reported gene set
job4_GSE131706_methyl.sh    450K array, SVZ
job5_GSE80017_methyl.sh     450K array, PFC
job6_GSE74193_methylarray.sh  450K array, schizophrenia
GSE*_run_analysis.R         minfi processing, quantile normalisation, limma on M-values
GSE131706_agecontrol.R      age as a covariate
GSE131706_sexcontrol.R      sex as a negative control
```

### Integration
```
chd8_asd_axis_v2.R              CHD8 targets against expression-matched nulls
chd8_dose_response_both_arms.R  nested set sizes, autism and schizophrenia
scz_methylation_axis.R          methylation-to-expression, directional subsets
ko_pathway_map.R                KEGG orthologue to pathway mapping
ko_pathway_DA_all_arms.R        microbial pathway differential abundance
kegg_axis_correlation.R         gut-brain pathway correlation, stratified
cross_arm_microbiome_test.R     cross-disorder genus overlap and SCFA enrichment
microbiome_epigenomics_axis.R   third axis: microbial against epigenomic
matched_methods_attempt.R       MOFA2 and DIABLO, documented failure
```

---

## Methods and tools

**Microbiome** — cutadapt · DADA2 · QIIME 2 · SILVA 138 (taxonomy) · SILVA 128 via SEPP (phylogeny) · MMUPHin · ANCOM-BC2 · MaAsLin2 · PICRUSt2 · MetaCyc · KEGG · UniFrac · PERMANOVA / PERMDISP

**Transcriptomics** — Trim Galore · STAR · featureCounts · DESeq2 · limma-voom · CAMERA · dtangle · ashr · apeglm

**Epigenomics** — minfi · limma · bowtie2 · MACS2 · csaw

Full citations for every tool are in the dissertation's reference list.

---

## Environment

`environment/` contains `sessionInfo()` output, the QIIME 2 version, the module
load lines used on the HPC, and `tool_versions.txt` — the exact module version of
every tool the job scripts loaded, read off the scripts themselves rather than
recorded by hand, annotated where one stage used a different version from the
rest.

Package versions differ between arms, and `tool_versions.txt` sets out how. The R
environment is the least uniform part of the project: three mechanisms were used,
not one. Most analyses ran against the module stack at R 4.4.1 with Bioconductor
3.19. The depression microbiome arm ran at R 4.3.1 with a user library ahead of
the module on the library path, which is why some of its packages resolve to
Bioconductor 3.20. And the differential-abundance steps ran under conda, outside
the module stack, because ANCOM-BC2 and MaAsLin2 would not resolve against it.

There are four conda environments across two project directories, not one shared
environment. The schizophrenia arm's ANCOM-BC2 ran in `microbiome_r`
(`ancombc2_scz.sh`); the depression arm's ANCOM-BC2, MaAsLin2 and pathway
differential abundance ran in a separately built `ancombc` environment; the
autism arm's leave-one-dataset-out check also used `microbiome_r`; and PICRUSt2
has one environment in each project.

This is a real limitation, and it falls on the result most exposed to it. The
differential-abundance tests for schizophrenia and depression were run under
independently built environments whose package versions were never pinned to each
other, so a version difference between them cannot be ruled out as a contributor
to a difference between those two arms — and schizophrenia and depression are
exactly the pair in which this project reports its one gut convergence. The
concordance requirement between ANCOM-BC2 and MaAsLin2 was applied within each
arm, not across them, so it does not address this. The convergence should be read
with that in mind.

`environment/env_setup/` holds the scripts that built the `ancombc` environment
and installed MMUPHin. The full account is in the dissertation's Supplementary
Methods S3.

---

## Known gaps

Two figures were produced in interactive sessions that were not saved: the
cross-arm microbiome figure and the cross-disorder set-level panel. Their inputs
are present and the values are reproducible from the saved tables, but the
plotting code is not in this repository.

The taxonomy-to-function bridge scripts are per-arm rather than a single
parameterised script; `bridge_asd.py` records in its header that it replicates the
schizophrenia logic, and the selection criteria are identical.

Two of the six manner-of-death IEG p-values have no script. `mdd_filtered_rerun.R`
computes GSE101521's suicide-against-control contrast and `mdd_filtered_complete.R`
the three GSE102556 contrasts, but nothing here computes the set-level p-values
for depression-without-suicide against controls, or for suicide against
non-suicide within depression. `deseq_suicide.R` runs both contrasts gene by gene
and stops there.

No script in this repository produces the corrected SCFA enrichment.
`RESULTS_MAP.md` cites `scfa_enrichment.csv` and `tested_genus_universe.csv` for
3 of 5 against 10 of 166, but `cross_arm_microbiome_test.R` writes only
`cross_arm_overlap.csv` and `shared_genera.csv`, and nothing here builds the
166-genus universe that replaced the earlier 231. That correction is recorded in
`RESULTS_MAP.md` but the code behind it is not here.

`GSE41826_run_analysis.R` is a placeholder, not a script. The file contains the
single line `[paste the entire script contents here]` and nothing else — it is 40
bytes — and the same is true of the copy it came from. GSE41826 is the depression
arm's primary methylation analysis, so the script behind that result is not here.
The dataset's download script and its sex control, `GSE41826_download.sh` and
`GSE41826_sexcontrol.R`, are present, and the sex control shows the data being
loaded; the primary analysis itself is not.

It appears to have been run interactively. All four HPC project directories were
searched for a script producing GSE41826's `mVals_betaVals.RData` and none
exists, and GSE41826
is the only one of the four array datasets with no `job*.sh` alongside `job4`,
`job5` and `job6`. The result stands on the saved objects and the dissertation's
methods; the code path to it is not recoverable from these directories.

Two further things the pipeline depends on are not in this repository. The first is
GSE74193's subject-resolution and prenatal-removal step, described below, which
has no script at all. The second is a script referenced by `RESULTS_MAP.md`:
`chd8_percohort.R`, cited for the per-cohort replication of the CHD8 consensus
set. Its log survives with the values in it, but the script has not been
recovered. Its input — the CHD8 peak-to-gene mapping — is produced by
`chd8_prededup.R`, which is present; the per-cohort step applied to that input is
not.

The two H3K9me2 window scripts are a sequence, not alternatives, and both are
required. `csaw_GSE215991.R` does the counting and writes `windows_raw.rds` and
`bins_10kb.rds`; `csaw_GSE215991_v2.R` reads those two files, does the testing,
and writes `GSE215991_csaw_v2.RData`, which the promoter-restricted test and both
controls then load. Run the first, then the second via `job_csaw_v2.sh`.

The two GSE74193 scripts are not interchangeable. They differ in where they get
their IDATs — `GSE74193_methylation_analysis.R` downloads them from GEO with
`getGEOSuppFiles`, `GSE74193_run_analysis.R` reads them from a local path — but
they also fit different models. `GSE74193_methylation_analysis.R` uses
`~0 + diagnosis + age + sex + period` with period as a factor.
`GSE74193_run_analysis.R` builds the design incrementally and adds period as
`as.numeric(period)`, which codes its three levels as 1, 2, 3 and so treats an
unordered categorical as a linear scale.

`job6_GSE74193_methylarray.sh` calls `run_analysis.R`, so which of the two ran is
on record — it is the one with period coded numerically.

The gap is elsewhere, and it applies whichever script ran: neither removes the
prenatal samples, and neither reduces the 675 arrays to one record per subject.
The 185 / 295 split in the dataset table depends on both steps, so both happened
outside these scripts — in something not in this repository, or interactively.
That filtering is the missing piece for this dataset, not the choice of script.

GSE131706 and GSE80017 come in the same pair, and split the same way: the
`_run_analysis.R` of each, called by `job4` and `job5`, builds its design
incrementally and leaves age out, while the `_methylation_analysis.R` of each
includes it. For these two the numeric coding is harmless — the variables added
that way, sex and region or tissue, have two levels each, so numeric and factor
coding give the same fit. `GSE131706_agecontrol.R` is the check that adds age
back. Only GSE74193's three-level period is affected by the coding, and only
GSE74193 has the filtering gap described above.

Several wrappers call an R script under a shorter name than the file carries
here, because the files were given accession prefixes when they were organised
for this repository. The wrapper is verbatim; the file is present under the
longer name.

| A wrapper calls | The repository holds |
|---|---|
| `agecontrol_ASD.R`, `sexcontrol_ASD.R` | `GSE131706_agecontrol.R`, `GSE131706_sexcontrol.R` |
| `sexcontrol_MDD.R` | `GSE41826_sexcontrol.R` |
| `promoter_test.R`, `sex_control.R`, `rizavi_genes.R` | `GSE215991_promoter_test.R`, `GSE215991_sex_control.R`, `GSE215991_rizavi_genes.R` |
| `run_analysis.R` | `GSE74193_run_analysis.R`, `GSE131706_run_analysis.R`, `GSE80017_run_analysis.R` — one per dataset, selected by the calling job script's path. GSE41826 has no job script and its `run_analysis.R` is the placeholder described above. |

featureCounts was not called identically across the arms. `--primary` is set in
the five schizophrenia counting scripts — `featurecounts_PRJNA314463`,
`featurecounts_PRJNA319583`, `featurecounts_PRJNA574470`,
`featurecounts_array_PRJNA574470` and `fc_heavy5` — and in
`featurecounts_GSE102741_perBAM`, and is absent everywhere else.

On the best available reading it made no difference to the counts. `-M` is not
set in any script here, and featureCounts counts multi-mapping reads only when
`-M` is given, so multi-mapping reads were excluded in every arm regardless of
`--primary`. The Subread documentation states that `-M` is ignored when
`--primary` is set, which if taken literally would mean the flag does change
behaviour; that line is widely reported as a documentation error. The
inconsistency is recorded because it is visible in the scripts. Whether it has
any effect has not been measured here, and counting one BAM both ways would
settle it.

`recover_PRJNA434002_6samples.sh` never ran, and is kept only as a record of an
attempt. It names six run accessions in the `SRR9262` series, but every library in
PRJNA434002 is `SRR9292` — as are the accessions named in this cohort's other
repair scripts, `star_2samples.sh` and `realign_and_fc.sh` — and the `SRR9262`
series appears nowhere else in the project. It also works in a different project
directory from the rest of the cohort's scripts. The job left no output: the STAR
results for this cohort predate it by a month. The cohort's count matrix holds 41
libraries, matching the deposited total, so nothing is missing from it. The script
also ends at STAR with `--quantMode GeneCounts` and never calls featureCounts,
which is consistent with its never having been completed.

Every other stage resolves to a script, with the four exceptions named above: GSE74193's subject resolution and prenatal removal, GSE41826's primary analysis, the corrected SCFA enrichment, and two of the six manner-of-death set-level p-values. PRJNA574470's counting looked ambiguous —
a serial and a job-array version of featureCounts with no way to choose between
them — until `fc_heavy5.sh` was recovered: the three run in sequence, the array
version superseding the serial one and `fc_heavy5.sh` counting the five heaviest
libraries separately before `merge_574470.sh` combines them. The detail is in
`02_transcriptomics/upstream/PROVENANCE.md`.

The job scripts are kept verbatim, with HPC paths, SLURM resource requests and
`module load` lines in place. They record what ran; they are not portable, and
the paths will need changing to run them anywhere else. Where a stage was
attempted more than once, the two `PROVENANCE.md` files say which attempt stands,
resolved against modification times rather than filenames — `_final` in a
filename turned out twice to mean the opposite.

---

## Citing the data

Deposits should be cited by their original publications, not by accession alone.
- **Autism microbiome** — PRJNA453621, PRJNA529598, PRJNA769228, PRJNA1022510
- **Schizophrenia microbiome** — PRJNA880407
- **Depression microbiome** — PRJNA591924, PRJNA776170
- **Autism transcriptomics** — GSE102741, PRJNA434002, GSE59288 / GSE51264 (Voineagu *et al.* 2011, *Nature* 474:380; Parikshak *et al.* 2016, *Nature* 540:423)
- **Schizophrenia transcriptomics** — GSE80655, GSE78936, GSE138082
- **Depression transcriptomics** — GSE101521, GSE102556, GSE80655
- **Autism epigenomics** — GSE61492 (Cotney *et al.* 2015, *Nat Commun* 6:6404), GSE131706, GSE80017 (Sun *et al.* 2022, *Psychiatr Genet* 32:131)
- **Schizophrenia epigenomics** — GSE74193, GSE215991
- **Depression epigenomics** — GSE41826 (Bustamante *et al.* 2018, *Psychiatry Res* 260:439)

Full citations for every deposit and every method are in the dissertation's
reference list. Where a deposit's source publication is not listed above, cite the
accession and the GEO or SRA record.

---

## Licence

MIT. See `LICENSE`.
