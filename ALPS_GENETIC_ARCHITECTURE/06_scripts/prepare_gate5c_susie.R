suppressPackageStartupMessages({
  library(data.table)
})

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5B")) && dir.exists(file.path(candidate, "Gate5C"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

root <- resolve_project_root()
harm_file <- file.path(root, "Gate5B/03_coloc_results/Gate5B_coloc_harmonized_regions.csv.gz")
bim_file <- Sys.getenv(
  "ALPS_LD_BIM",
  unset = file.path(dirname(root), "04_ld_reference", "1000G.EUR.QC.bim")
)
out_dir <- file.path(root, "Gate5C/05_finemapping")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(harm_file)) stop("Missing harmonized regional input: ", harm_file)
if (!file.exists(bim_file)) stop("Missing EUR LD BIM file; set ALPS_LD_BIM: ", bim_file)

harm <- fread(cmd = paste("gzip -dc", shQuote(harm_file)))
harm <- harm[
  region_id == "ALPS_L16_011" &
  external_trait == "WMH" &
  ALPS_phenotype %in% c("mALPS", "left_ALPS", "right_ALPS")
]
harm_valid <- harm[
  is.finite(beta_alps) & is.finite(se_alps) &
  is.finite(beta_ext_aligned) & is.finite(se_ext)
]
bim <- fread(
  bim_file,
  col.names = c("CHR", "SNP", "BP_ref", "ALT_ref", "REF_ref"),
  select = c(1, 2, 4, 5, 6)
)
bim <- bim[CHR == 16 & BP_ref >= 86736383 & BP_ref <= 87736383]
bim <- bim[!duplicated(SNP)]

common <- Reduce(
  intersect,
  c(
    lapply(split(harm_valid$SNP, harm_valid$ALPS_phenotype), unique),
    list(bim$SNP)
  )
)
common_bim <- bim[SNP %in% common][order(BP_ref)]
if (nrow(common_bim) < 100L) stop("Too few common SNPs for SuSiE preparation: ", nrow(common_bim))
harm_common <- merge(
  harm_valid[SNP %in% common],
  common_bim,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)
if (anyDuplicated(common_bim$SNP)) stop("Duplicate SNPs remain in the LD reference subset.")
if (anyNA(harm_common$REF_ref)) stop("Reference alleles are missing after BIM alignment.")

fwrite(common_bim[, .(SNP)], file.path(out_dir, "susie_common_snps.txt"), col.names = FALSE)
fwrite(harm_common, file.path(out_dir, "susie_harmonized_common.csv.gz"))

audit <- data.table(
  metric = c(
    "common_SNPs_all_three_ALPS_WMH_reference",
    "reference_region_SNPs",
    "mALPS_pair_SNPs_before_reference",
    "left_pair_SNPs_before_reference",
    "right_pair_SNPs_before_reference"
  ),
  value = c(
    nrow(common_bim),
    nrow(bim),
    uniqueN(harm[ALPS_phenotype == "mALPS", SNP]),
    uniqueN(harm[ALPS_phenotype == "left_ALPS", SNP]),
    uniqueN(harm[ALPS_phenotype == "right_ALPS", SNP])
  )
)
fwrite(audit, file.path(out_dir, "susie_input_audit.csv"))
print(audit)
