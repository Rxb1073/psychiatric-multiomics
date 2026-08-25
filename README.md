# Multi-omic reanalysis of autism, schizophrenia and major depression

Analysis code for an MSc Bioinformatics dissertation. Twenty-three open-access
deposits were retrieved and twenty-one analysed across three molecular layers —
gut microbial 16S rRNA amplicon sequencing, post-mortem cortical transcriptomics,
and post-mortem cortical DNA methylation and histone modification — with each
disorder processed through the same pipeline so that differences between them
reflect biology rather than the analytical choices of separate published studies.

**No data are included in this repository.** All deposits are public and listed
below with their accessions. Sample sheets are excluded because they carry
donor-level covariates.

---

## What was found

Three disorder-specific brain signatures, none extending beyond its own disorder:
GABAergic interneuron markers reduced in schizophrenia, immediate-early genes in
depression, and CHD8-bound promoters displaced in autism. One convergence, in the
gut, between schizophrenia and depression.

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
| | GSE78936 | BA11 | poly(A), PE | 52 samples | 40 donors (BA11 28) | 16 / 12 | primary (frontal) | −s 0; 12 donors at BA24 and BA9; BA11 subset never duplicated |
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
```
qiime_phase1.sh          import, primer removal, DADA2 denoising per cohort
qiime_phase2.sh          taxonomy (SILVA 138), genus collapse, prevalence filter
qiime_sepp.sh            fragment insertion into the SILVA 128 reference tree
stage3_mmuphin.R         batch correction, diagnosis retained as covariate
stage4_diversity.R       alpha and beta diversity, PERMANOVA with PERMDISP
stage5_diffabund.R       ANCOM-BC2 and MaAsLin2; concordance required
lodo_check.R             leave-one-dataset-out
picrust2_run.sh          predicted gene content, MetaCyc and KEGG
bridge_*.py              taxonomy-to-function bridge, one per arm
```

### Transcriptomics
```
trim_*.sh, fc_*.sh       trimming, STAR alignment, featureCounts
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
minfi processing, quantile normalisation, limma on M-values
CHD8 peak calling (MACS2), three-antibody consensus, CHD7 specificity control
H3K9me2 csaw sliding windows, genome-wide and promoter-restricted
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

`environment/` contains `sessionInfo()` output, the QIIME 2 version, and the
module load lines used on the HPC. Package versions differ between arms; those
departures are documented in the dissertation's Supplementary Methods S3 and
summarised in `environment/departures.md`.

---

## Known gaps

Two figures were produced in interactive sessions that were not saved: the
cross-arm microbiome figure and the cross-disorder set-level panel. Their inputs
are present and the values are reproducible from the saved tables, but the
plotting code is not in this repository.

The taxonomy-to-function bridge scripts are per-arm rather than a single
parameterised script; `bridge_asd.py` records in its header that it replicates the
schizophrenia logic, and the selection criteria are identical.

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
