# DTI-ALPS–WMH reproducibility release

This directory is the publication-ready index for the analysis supporting:

> Genetic architecture of DTI-ALPS phenotypes and shared genetic links with cerebral small vessel disease

It does not redistribute third-party GWAS summary statistics. The data manifest provides source URLs, identifiers, phenotype definitions, sample-overlap notes, and the local files used in the analysis. Users must obtain source data under the terms of the originating studies.

## Analysis order

1. `06_scripts/build_alps_architecture.py` — ALPS input QC, significant variants, clumping, and merged loci.
2. `03_LDSC/` and `06_scripts/parse_ldsc_results.py` — ALPS SNP heritability.
3. `06_scripts/compile_gate5a_phewas.py` — targeted cross-trait screen.
4. `06_scripts/run_gate5c_ldsc_rg.sh` and `run_gate5c_ldsc_rg_lateral.sh` — cross-trait LDSC.
5. `06_scripts/run_gate5b_coloc.R` — classic regional colocalization.
6. `06_scripts/prepare_gate5c_susie.R` and `run_gate5c_susie_coloc.R` — LD alignment, SuSiE, and multi-signal colocalization.
7. `06_scripts/run_gate5c_independent_wmh_coloc.R` and `audit_charge_eur_wmh_replication.py` — non-UK Biobank regional evaluations.
8. `06_scripts/extract_gate5d_credible_variants.R`, `query_gate5d_gtex_vep.py`, and `run_gate5d_magma.sh` — functional annotation.
9. `06_scripts/run_malps_to_nonukb_wmh_mr.R` — exploratory supplementary MR.
10. `09_manuscript/scripts/build_main_figures_tables.R` — main tables and figures.
11. `09_manuscript/scripts/build_figure5_q2.R` — redesigned data-selected chr16 regional figure and full coloc-prior sensitivity.
12. `09_manuscript/scripts/build_embedded_word_document.py` — manuscript Word assembly.

All paths above are relative to the `ALPS_GENETIC_ARCHITECTURE` project root.

## Locked primary parameters

- Genome-wide significance: `P < 5 × 10^-8`
- LD clumping: 1000 Genomes European reference, `r² < 0.001`, 10 Mb
- Cross-trait locus interval: ±500 kb
- Data-selected chr16 interval: GRCh37 `chr16:86,736,383–87,736,383`
- Classic coloc: `p1 = 1 × 10^-4`, `p2 = 1 × 10^-4`, primary `p12 = 1 × 10^-5`
- Coloc sensitivity: `p12 = 1 × 10^-4, 1 × 10^-5, 1 × 10^-6, 1 × 10^-7`
- Strong model-based regional support: `PP4 > 0.80`
- SuSiE: `L = 10`, 95% credible-set coverage, maximum 500 iterations
- LDSC genetic-correlation multiplicity: Bonferroni `P < 0.00278` and BH-FDR `< 0.05`

## Interpretation locks

- chr16 was selected after the primary WMH screen. Regional analyses are susceptible to winner selection.
- The primary ALPS and WMH GWAS substantially overlap. The targeted screen, LDSC, coloc, and SuSiE are not independent replications.
- SuSiE handles multiple association signals but does not correct cross-trait sampling covariance from participant overlap.
- Formal Traylor colocalization did not support sharing; CHARGE could not support posterior estimation.
- DTI-ALPS is an orientation-dependent periventricular diffusion phenotype, not a direct measurement of glymphatic clearance.
- MR is exploratory complementary directional evidence and is not external replication or proof of causality.

## Release procedure

Before submission:

1. Run the manuscript and figure QA scripts.
2. Run `generate_checksums.py` from this directory.
3. Remove private absolute paths from release copies of logs.
4. Publish the curated code, parameter locks, non-restricted derived tables, and figure source data to GitHub or OSF.
5. Archive the exact release on Zenodo and insert the repository URL and DOI in the manuscript.

The repository URL and DOI remain an author-controlled publication action and must not be fabricated.
