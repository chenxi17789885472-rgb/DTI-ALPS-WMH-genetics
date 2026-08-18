# R code reproducibility audit

Audit date: 2026-08-13

## Scope

Core R analyses: conventional colocalization, non-UK Biobank colocalization, SuSiE multi-signal fine-mapping/colocalization, credible-set extraction, exploratory MR, and final figure/table generation.

## Corrections made before release

1. Replaced hard-coded local project paths with automatic project-root discovery and optional `ALPS_PROJECT_ROOT` override.
2. Added `ALPS_LD_BIM` override and explicit checks for the external EUR LD BIM file.
3. Added required-column, finite-value, positive-standard-error, duplicate-SNP, instrument-strength, and minimum-SNP checks.
4. Added LD dimension, diagonal, symmetry, finiteness, and range checks.
5. Fixed and documented random seeds for conventional coloc bookkeeping, SuSiE, and MR-PRESSO.
6. Added explicit SuSiE convergence gates and a software-version log.
7. Removed the Figure 1 script's ad hoc local-library precedence so that it uses the frozen renv environment.
8. Excluded the obsolete R audit script tied to the superseded seven-figure/three-table layout.

## Execution results

- R syntax checks: PASS for all included R scripts.
- Conventional coloc rerun: PASS.
- Traylor non-UK Biobank coloc rerun: PASS.
- SuSiE rerun: PASS; WMH, mean ALPS, left ALPS, and right ALPS models all converged.
- Credible-set extraction rerun: PASS.
- Exploratory MR and sensitivity workflow rerun: PASS.
- Final figure/table rebuild: PASS.
- Standalone package execution with `run_core_reproduction.R --execute`: PASS.
- Figure/table/supplement reconciliation: 93/93 checks passed.
- `renv::status()`: No issues found.

## Exact key-output reproducibility

The following SHA-256 hashes were identical before and after code correction and rerun:

- `Gate5B_coloc_summary.csv`: `faf1a1637e35647913ea99e6a3d51e88248427e9c476f6e4b96e1fc3004d75e9`
- `Gate5C_independent_WMH_coloc_summary.csv`: `d198541255b0c5f01be9c2740833a758eb11577b8a7ef7c3c7d8fc713ac721c2`
- `Gate5C_SuSiE_coloc_signal_pairs.csv`: `eadd2b1ad73b04c21eca28ec009fb96443fb80274791e7036977beeb40d47395`
- `Gate5C_SuSiE_credible_sets.csv`: `2149c6ca24c689ac2769f021a104b4e2553923e207b5603f2713a43321c79db1`
- `mALPS_to_Traylor2016_WMH_MR_results.csv`: `412e490dc3974d2efc68157925d064262b7ed6a855e5da3ec43f312431964cb9`

## Remaining boundary

This audit establishes computational reproducibility for the supplied inputs and locked environment. It does not convert overlap-sensitive regional evidence into independent validation, nor does it resolve biological causality. Full LDSC and MAGMA regeneration also requires their external reference files and full GWAS inputs, which are indexed but not redistributed here.
