# DTI-ALPS–WMH genetic relationship: reproducibility materials

This repository accompanies the manuscript **"Genome-wide and Regional Genetic Relationships Between DTI-ALPS Phenotypes and White Matter Hyperintensity Burden: Colocalization and Fine-Mapping Analyses."** It reproduces the core R-based analyses and publication figures. It contains no complete third-party GWAS dataset; public source data must be obtained from the original repositories under their respective terms.

## Authors

Qingqing Fu and Ge Zhang contributed equally. Weiyuan Huang and Feng Chen are
co-corresponding authors. All authors are affiliated with the Department of
Radiology, Hainan General Hospital, Affiliated Hospital of Hainan Medical
University, Haikou, Hainan, China.

## Repository scope

The release contains analysis scripts, frozen software metadata, regional
analysis inputs required by the executable workflows, derived numerical
outputs, figure source data, and publication figures. It does not contain
individual-level participant data. See `LICENSE.md` and
`ALPS_GENETIC_ARCHITECTURE/09_manuscript/reproducibility_release/DATA_SOURCES.tsv`
for the licensing and provenance boundaries.

## Verified environment

- R 4.6.0 (aarch64-apple-darwin23)
- Package versions are frozen in `renv.lock`.
- At release, `renv::status()` returned `No issues found`.

Restore the R environment from this directory:

```r
install.packages("renv")
renv::restore()
```

## Quick validation

From the package root:

```bash
Rscript run_core_reproduction.R
```

This checks the R scripts, required files, frozen packages, and locked analysis parameters without overwriting results.

Release-file integrity can be checked against `SHA256SUMS.txt`. Regenerate the
manifest after an intentional release update with:

```bash
python3 generate_release_checksums.py
```

To rerun the core R analyses and publication figures:

```bash
Rscript run_core_reproduction.R --execute
```

The executable sequence is:

1. conventional regional colocalization;
2. non-UK Biobank Traylor colocalization;
3. SuSiE multi-signal fine-mapping and colocalization using the frozen regional LD matrix;
4. credible-set extraction;
5. exploratory mALPS-to-WMH MR;
6. Figure 1 and Figures 2–6, Tables 1–2, and Supplementary Figure S1.

## Important data boundary

The package includes the extracted chromosome 16 Traylor file, harmonized regional inputs, the frozen 1,832-variant LD matrix, derived tables, and figure source data. It does not include the 193-MB genome-wide 1000 Genomes BIM file or full third-party GWAS summary statistics. Therefore, `prepare_gate5c_susie.R` requires either the original workspace layout or the `ALPS_LD_BIM` environment variable. The supplied frozen regional inputs are sufficient to rerun `run_gate5c_susie_coloc.R` directly.

## Interpretation locks

- Primary ALPS and WMH datasets overlap substantially; primary regional analyses are internal, overlap-sensitive evidence.
- The chromosome 16 interval was data-selected and is subject to winner selection.
- Traylor formal colocalization did not support sharing; CHARGE was not suitable for posterior estimation.
- DTI-ALPS is an MRI-derived diffusion phenotype, not a direct measure of glymphatic clearance.
- MR is complementary exploratory directional evidence, not proof of causality or external replication.

## Code audit status

All included R scripts passed syntax checks. Core colocalization, SuSiE, credible-set, MR, and figure workflows were rerun after path-portability and input-validation changes. Five key result files retained identical SHA-256 hashes. See `R_CODE_REPRODUCIBILITY_AUDIT.md`.

The obsolete `audit_submission_assets.R`, which expected the superseded seven-figure/three-table layout, is intentionally excluded. The current manuscript uses Figures 1–6, Tables 1–2, Supplementary Figure S1, and Supplementary Tables S1–S13.

## Citation and permanent archive

Citation metadata are provided in `CITATION.cff`. The version-specific Zenodo
DOI will be added after the public GitHub release has been archived.

## License

Original source code is available under the MIT License. Third-party data,
database responses, and externally sourced images retain their original terms.
See `LICENSE.md` for details.
