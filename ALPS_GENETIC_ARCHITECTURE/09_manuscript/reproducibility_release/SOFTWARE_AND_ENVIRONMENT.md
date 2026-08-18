# Software and environment

The final manuscript build was prepared on macOS in July 2026.

Core analysis software:

- R 4.6.0
- coloc 5.2.3
- susieR 0.14.2
- TwoSampleMR 0.7.9
- MR-PRESSO 1.0
- MAGMA 1.10
- PLINK 2.0
- LDSC (project-pinned source under `03_LDSC/software/ldsc-ldsc39`)

Figure and document software:

- ggplot2 4.0.3
- patchwork 1.3.2
- data.table
- ggrepel
- openxlsx
- python-docx
- openpyxl

The exact source files and analysis-specific version records are retained in:

- `Gate5B/05_logs/software_versions.txt`
- `Gate5C/07_logs/independent_WMH_coloc_software_versions.txt`
- `Gate5C/10_mALPS_to_WMH_MR/05_logs/software_versions.txt`
- `Gate5C/05_finemapping/Gate5C_SuSiE_audit.csv`

Before public archiving, run `sessionInfo()` and package-version capture in the final environment and append the output to this file.
