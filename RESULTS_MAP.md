# Results map

Every quantitative claim in the dissertation, against the script that produced it
and the output file that holds it. Paths are relative to the project root on the
HPC; the scripts are in this repository under the corresponding layer folder.

This file exists because a number that cannot be traced to a file is a number that
cannot be defended. Several figures were corrected during writing when this
mapping was built — those corrections are noted where relevant.

---

## Transcriptomics · schizophrenia

| Claim | Value | Script | Output |
|---|---|---|---|
| Frontal donors after resolution | 109 samples → 76 donors | `scz_pooled_frontal.R` | `Transcriptomics_SCZ/tables/` |
| GABAergic set-level, unadjusted | FDR 8.2 × 10⁻³, ρ = 0.389 | `scz_camera_igc_check.R` | `CELLTYPE_SETS/camera_igc/CAMERA_estimated.csv` |
| GABAergic, composition-adjusted | FDR 1.0 × 10⁻³ (MB neuronal), 1.4 × 10⁻³ (CA neuronal), 1.4 × 10⁻³ (CA inhibitory) | `scz_camera_composition_adjusted.R` | `Transcriptomics_SCZ/tables/camera_composition/CAMERA_composition_comparison.csv` |
| Within-set correlation under adjustment | 0.389 → 0.100 | same | same |
| Immediate-early, unadjusted / adjusted | FDR 0.112 / 0.034 | same | same |
| Inhibitory fraction reduced | P = 0.008 (Wilcoxon) | `scz_composition_alltypes.R` | `Transcriptomics_SCZ/tables/SCZ_composition_tests_alltypes.csv` |
| Total neuronal fraction reduced | P = 0.040 | same | same |
| Microglial proportion unchanged | P = 0.914 (Cell Atlas), 0.645 (MultiBrain) | same | same |
| Source-publication re-analysis | 0, 7 and 0 genes | — | recorded in `Transcriptomics_SCZ/README.txt` |

**Correction made during verification.** The chapter's inhibitory and neuronal
P values (0.008, 0.040) come from an unadjusted Wilcoxon test, not from the linear
model adjusting for dataset. `scz_composition_alltypes.R` reports both so the two
are not conflated. No result in the report mixes them.

---

## Transcriptomics · depression

| Claim | Value | Script | Output |
|---|---|---|---|
| Donors, primary pool | 167 → 107 | `mdd_primary_2cohort.R` | `Transcriptomics_MDD/tables/01_PRIMARY_pooled_DE.tsv` |
| Immediate-early module | FDR 0.023 | `mdd_filtered_complete.R` | `Transcriptomics_MDD/tables/` |
| ERK module | FDR 0.023 | same | same |
| Robustness, weakest of four | P = 2.3 × 10⁻² | same | `10_robust_leave_one_gene_out.tsv`, `11_robust_leave_one_donor_out.tsv` |
| Within-cohort permutation | 5,000 label shuffles | `mdd_ieg_final_checks.R` | logged |
| Manner-of-death contrasts | six rows, see below | three GSE102556 contrasts `mdd_filtered_complete.R`; GSE101521 suicide vs control `mdd_filtered_rerun.R` §6 (`grpMDDs - grpCON`). The other two GSE101521 IEG p-values have no script — see Known gaps | `16_suicide_decomposition.tsv` (`suicide/suicide_decomposition.tsv` as written) holds the GSE102556 rows only |

**Manner-of-death subgroup sizes**, verified against
`sheets/PRJNA398031_SraRunTable.csv` (GSE102556) and the `group3` column of
`sheets/GSE101521_colData.tsv`:

| Contrast | n | IEG P |
|---|---|---|
| GSE102556 · suicide vs non-suicide, within controls | 11 v 11 | 0.296 |
| GSE102556 · MDD vs suicide controls | 26 v 11 | 0.031 |
| GSE102556 · MDD vs non-suicide controls | 26 v 11 | 0.692 |
| GSE101521 · depression without suicide | 9 v 29 | 0.258 |
| GSE101521 · depression with suicide | 21 v 29 | 0.0040 |
| GSE101521 · suicide vs non-suicide, within depression | 21 v 9 | 0.449 |

All 26 GSE102556 MDD donors died by suicide (143 of 143 deposit-wide), so the
second row holds manner of death constant. The 11 non-suicide controls comprise
7 natural and 4 accident deaths.

---

## Transcriptomics · autism

| Claim | Value | Script | Output |
|---|---|---|---|
| Donors after age matching | 45 (16 ASD / 29 control) | `asd_child_confoundfree.R` | `Transcriptomics_ASD/tables/` |
| Genes at FDR | 1 (*CXCL8*), 0 under random effects | `asd_metaanalysis.R` | `DE_ASD_child_metaanalysis_confoundfree.csv` |
| Marker sets, no composition covariates | IEG FDR 0.219, GABAergic 0.179 | `asd_setlevel_camera.R` | `Transcriptomics_ASD/tables/setlevel/ASD_setlevel_summary.csv` |
| Marker sets, with composition | IEG 0.373, GABAergic 0.382 | same | same |
| Nominal excess over chance | 1,780 of 17,485 against ~874 | `asd_metaanalysis.R` | logged |
| Autism composition tested separately | Astro β = 0.043, P = 0.008; Neuro β = −0.052, P = 0.101 | `02_transcriptomics/asd/` (post-hoc, 45 donors) | `Transcriptomics_ASD/tables/nocomposition/ASD_composition_tested_separately.csv` |

---

## Epigenomics

| Claim | Value | Source |
|---|---|---|
| Schizophrenia arrays → subjects | 675 → 515 → 480 postnatal | `Transcriptomics_SCZ` methylation scripts |
| Differentially methylated positions | 71,877 at FDR; 213 after \|Δβ\| ≥ 0.05 | `DMCpGs_SCZ_postnatal_autosomal_213.csv` |
| Calibration ratios | 94.0, 0.54, 0.79, 0.93, 0.33 | Table 4 of the dissertation |
| H3K9me2 | 0 of 30,792 genome-wide; 0 of 49,904 promoter | csaw output |
| CHD8 consensus sites | 6,015 shared → 508 CHD7-matched removed → 5,507 | peak-calling scripts |

---

## Integration

| Claim | Value | Script | Output |
|---|---|---|---|
| CHD8 sites within 2 kb of a TSS | 4,601 of 5,507 (83.5%) | `chd8_prededup.R` | `07_chd8_axis/CHD8_peak_to_gene_prededup.csv` |
| Unique target genes | 4,351 | `chd8_asd_axis_v2.R` | `07_chd8_axis/CHD8_peak_to_gene.csv` |
| Targets with a meta-analytic estimate | 3,820 | same | same |
| Abundance ratio, targets vs non-targets | 2.8–3.6× | same | logged |
| Displacement at four set sizes | +0.0277, +0.0253, +0.0218, +0.0194 | `chd8_dose_response_both_arms.R` | `13_epigenomic_cross_disorder/CHD8_dose_response_both_arms.csv` |
| Matched nulls | −0.0019, −0.0040, −0.0055, −0.0073 | same | same |
| Schizophrenia comparison | six- to eightfold smaller; P = 0.434, 0.286, 0.086, 0.045 | same | same |
| Consensus set, per-cohort replication | +0.0283 and +0.0172, P = 1 × 10⁻⁴ | `chd8_percohort.R` | `07_chd8_axis/chd8_percohort.log`, `CHD8_percohort_matched.csv` |
| Sign test | 57.2% up, P < 2 × 10⁻¹⁶ | same | same |
| Heterogeneity | mean I² 9.8 vs 13.5 | same | same |
| Methylation subsets | promoter-hypo n = 44, P = 0.038; body-hypo P = 0.977; undirected P = 0.025 | `scz_methylation_axis.R` | `08_methylation_axis/` |
| Gut–brain correlation | ρ = −0.214 (P = 0.019), −0.228 (P = 0.017) | `kegg_axis_correlation.R` | `06_kegg_axis/cross_omics/` |
| Stratified | −0.152 and +0.049 | same | same |
| Distribution asymmetry | 78.2% vs 45.4% negative | same | same |
| Third axis, autism | 121 testable, 24 and 6 significant, 0 overlap, P = 1.0 | `microbiome_epigenomics_axis.R` | `14_microbiome_epigenomics/` |
| MOFA2 / DIABLO failure | 0 shared identifiers, 76 v 100 | `matched_methods_attempt.R` | logged |

**Correction made during verification.** The null medians were originally taken
from a 2,000-permutation run while the p-values came from a 10,000-permutation
run. All values now come from the 10,000 run, and Supplementary Table S2 reports
them consistently.

---

## Microbiome

| Claim | Value | Script | Output |
|---|---|---|---|
| Concordant genera | ASD 16 taxa (12 named), SCZ 26, MDD 4 | `stage5_diffabund.R`, `ancombc2_scz.R`, `ancombc_mdd_conda.R` | `stage5_concordant_core.csv`, `robust_genera_both_methods_v2.tsv`, `robust_genera_both_methods.tsv` |
| SCZ–MDD shared genera | 3 of a pool of 116, P = 0.035, all three same direction | `cross_arm_microbiome_test.R` | `12_cross_arm_microbiome/cross_arm_overlap.csv` |
| SCFA enrichment, five pairwise-shared | 3 of 5 against 10 of 166, 9.96-fold, P = 0.0015 | same | `12_cross_arm_microbiome/scfa_enrichment.csv` |
| Tested genus universe | 166 | same | `12_cross_arm_microbiome/tested_genus_universe.csv` |
| *Turicibacter* | +0.307 (ASD, MaAsLin2 q = 1.6 × 10⁻³); −1.207 (SCZ, ANCOM-BC2 q = 4.4 × 10⁻⁷) | — | `CROSS_ARM_figures/FigX_source_data.csv` |
| Fermentation pathways, SCZ | CENTFERM q = 0.0042, PWY-6590 q = 0.0042, FERMENTATION q = 0.0433 | `pathway_DA` | `Microbiome_SCZ/tables/pathway_SZ_vs_HC.tsv` |
| SCFA pathways, ASD and MDD | none significant; closest q = 0.18 and 0.34 | same | `Microbiome_ASD/tables/`, `Microbiome_MDD/tables/` |
| NSTI | 0.029–0.080 across cohorts | PICRUSt2 | `*/picrust2_out/KO_metagenome_out/weighted_nsti.tsv.gz` |
| Species-level resolution | 12.3% of 2,490 ASVs carry a binomial | — | taxonomy tables |

**Correction made during verification.** The SCFA background was originally
computed against a 231-genus universe that corresponded to no saved file. The
correct universe is the 166 genera in `tested_genus_universe.csv`. A restricted
schizophrenia–depression enrichment was also computed but not retained, because
its pool size depended on which method's tested set was used and no version was
written to file.

---

## Known gaps

The cross-arm microbiome figure and the cross-disorder set-level panel were drawn
in interactive sessions that were not saved. Their input values are in
`CROSS_ARM_figures/FigX_source_data.csv` and the CAMERA output tables
respectively, so the figures are reproducible from saved data, but the plotting
code is not in this repository.

---

## Microbiome · beta diversity

Recomputed after the first submission so that all three arms report the same
statistics from the same function.

| Claim | Value | Script | Output |
|---|---|---|---|
| Schizophrenia, five metrics | R² 0.047–0.078, all P = 0.001 | `01_microbiome/scz_beta_diversity_adonis2.R` | `Microbiome_SCZ/tables/SCZ_beta_diversity_adonis2.csv` |
| Depression, four metrics | R² 0.007–0.011, P = 0.005–0.099 | `01_microbiome/mdd_beta_diversity_adonis2.R` | `Microbiome_MDD/tables/MDD_beta_diversity_adonis2.csv` |
| PERMDISP, both arms | reported alongside each PERMANOVA | same scripts | same files |

Depression models include study as a covariate and report the marginal effect of
diagnosis, so its R² is not directly comparable with schizophrenia's single-cohort
value. Autism and depression Aitchison and Bray–Curtis values come from the
original `stage4_diversity.R` and `diversity_mdd.R`.


---

## A note on file names

Output paths in this table are the names the files carry in the dissertation's
results directory, not always the names the scripts write. The depression arm's
tables were copied and renumbered into that directory by a packaging step
(`build_mdd_v2.sh`, not included here, since it moves files rather than analysing
them). So `01_PRIMARY_pooled_DE.tsv` is `pooled_2cohort.tsv` as written by
`mdd_primary_2cohort.R`, `10_robust_leave_one_gene_out.tsv` is
`robust/leave_one_gene_out.tsv`, and so on. Where this table names a numbered
file, the script that produced it is the one in the Script column under its own
output name.

`robust_genera_both_methods.tsv` and `robust_genera_both_methods_v2.tsv` are
written by `concordance.R` and `concordance_v2.R`, not by the scripts cited in
the rows that consume them.

No script in this repository writes a separate manner-of-death table. Of the six
contrasts, four have a script: the three GSE102556 contrasts in
`mdd_filtered_complete.R` and GSE101521's suicide-against-control in
`mdd_filtered_rerun.R` §6. `deseq_suicide.R` runs all three GSE101521 contrasts
gene by gene but computes no set-level p-value, so the IEG figures for depression
without suicide against controls (9 v 29) and for suicide against non-suicide
within depression (21 v 9) are not produced by anything here.
