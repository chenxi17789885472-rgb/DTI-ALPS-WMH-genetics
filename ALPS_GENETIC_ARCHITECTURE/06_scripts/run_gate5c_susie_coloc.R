suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})

set.seed(20260725)

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5C", "05_finemapping"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

root <- resolve_project_root()
fm <- file.path(root, "Gate5C/05_finemapping")
harm_file <- file.path(fm, "susie_harmonized_common.csv.gz")
vars_file <- file.path(fm, "1000G_EUR_chr16_susie.unphased.vcor1.vars")
ld_file <- file.path(fm, "1000G_EUR_chr16_susie.unphased.vcor1")

snps <- fread(vars_file, header = FALSE)[[1]]
ld <- as.matrix(fread(ld_file, header = FALSE))
storage.mode(ld) <- "double"
rownames(ld) <- colnames(ld) <- snps
if (anyDuplicated(snps)) stop("Duplicate SNP identifiers in the LD variant list.")
if (nrow(ld) != ncol(ld) || nrow(ld) != length(snps) ||
    any(!is.finite(ld)) || max(abs(diag(ld) - 1)) > 1e-8 ||
    max(abs(ld - t(ld))) > 1e-8 || max(abs(ld)) > 1 + 1e-8) {
  stop("Invalid LD matrix dimensions or diagonal.")
}

harm <- fread(cmd = paste("gzip -dc", shQuote(harm_file)))
harm[, `:=`(
  EA_alps = toupper(EA_alps),
  OA_alps = toupper(OA_alps),
  REF_ref = toupper(REF_ref)
)]

orient_to_ref <- function(dt, beta_col) {
  b <- dt[[beta_col]]
  out <- fifelse(
    dt$EA_alps == dt$REF_ref,
    b,
    fifelse(dt$OA_alps == dt$REF_ref, -b, NA_real_)
  )
  out
}

make_dataset <- function(pheno, external = FALSE) {
  dt <- harm[ALPS_phenotype == pheno]
  dt <- dt[match(snps, SNP)]
  if (anyNA(dt$SNP)) stop("Missing summary statistics after LD ordering: ", pheno)
  if (external) {
    beta <- orient_to_ref(dt, "beta_ext_aligned")
    se <- dt$se_ext
    n <- 48454
  } else {
    beta <- orient_to_ref(dt, "beta_alps")
    se <- dt$se_alps
    n <- 31021
  }
  if (anyNA(beta) || anyNA(se)) stop("Allele mismatch or missing beta/SE: ", pheno)
  list(
    beta = beta,
    varbeta = se^2,
    snp = snps,
    type = "quant",
    N = n,
    sdY = 1,
    LD = ld
  )
}

extract_cs <- function(s, trait) {
  cs_list <- s$sets$cs
  if (is.null(cs_list) || length(cs_list) == 0) {
    return(data.table(
      trait = trait, signal = NA_character_, lead_SNP = NA_character_,
      lead_PIP = NA_real_, credible_set_size = 0L,
      coverage = NA_real_, min_abs_corr = NA_real_,
      mean_abs_corr = NA_real_, median_abs_corr = NA_real_
    ))
  }
  purity <- as.data.table(s$sets$purity)
  rbindlist(lapply(seq_along(cs_list), function(i) {
    idx <- cs_list[[i]]
    signal_name <- names(cs_list)[i]
    component <- as.integer(sub("^L", "", signal_name))
    lead_idx <- which.max(s$alpha[component, ])
    data.table(
      trait = trait,
      signal = signal_name,
      lead_SNP = snps[lead_idx],
      lead_PIP = s$pip[lead_idx],
      credible_set_size = length(idx),
      coverage = s$sets$coverage[i],
      min_abs_corr = purity[i, min.abs.corr],
      mean_abs_corr = purity[i, mean.abs.corr],
      median_abs_corr = purity[i, median.abs.corr]
    )
  }))
}

message("Running SuSiE for WMH...")
d_wmh <- make_dataset("mALPS", external = TRUE)
check_dataset(d_wmh, req = "LD")
s_wmh <- runsusie(
  d_wmh, L = 10, coverage = 0.95, maxit = 500,
  estimate_prior_variance = TRUE
)
if (!isTRUE(s_wmh$converged)) stop("WMH SuSiE model did not converge.")
saveRDS(s_wmh, file.path(fm, "susie_WMH.rds"))

all_cs <- list(extract_cs(s_wmh, "WMH"))
all_coloc <- list()
all_pip <- list()

for (pheno in c("mALPS", "left_ALPS", "right_ALPS")) {
  message("Running SuSiE for ", pheno, "...")
  d_alps <- make_dataset(pheno, external = FALSE)
  check_dataset(d_alps, req = "LD")
  s_alps <- runsusie(
    d_alps, L = 10, coverage = 0.95, maxit = 500,
    estimate_prior_variance = TRUE
  )
  if (!isTRUE(s_alps$converged)) stop(pheno, " SuSiE model did not converge.")
  saveRDS(s_alps, file.path(fm, paste0("susie_", pheno, ".rds")))
  all_cs[[length(all_cs) + 1]] <- extract_cs(s_alps, pheno)

  coloc_res <- coloc.susie(
    s_alps, s_wmh,
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
  )
  coloc_tab <- as.data.table(coloc_res$summary)
  coloc_tab[, `:=`(
    ALPS_phenotype = pheno,
    external_trait = "WMH",
    max_PP4_for_ALPS = max(PP.H4.abf, na.rm = TRUE)
  )]
  all_coloc[[length(all_coloc) + 1]] <- coloc_tab

  all_pip[[length(all_pip) + 1]] <- data.table(
    ALPS_phenotype = pheno,
    SNP = snps,
    PIP_ALPS = s_alps$pip,
    PIP_WMH = s_wmh$pip
  )
}

cs_tab <- rbindlist(all_cs, fill = TRUE)
coloc_tab <- rbindlist(all_coloc, fill = TRUE)
pip_tab <- rbindlist(all_pip, fill = TRUE)

fwrite(cs_tab, file.path(fm, "Gate5C_SuSiE_credible_sets.csv"))
fwrite(coloc_tab, file.path(fm, "Gate5C_SuSiE_coloc_signal_pairs.csv"))
fwrite(pip_tab, file.path(fm, "Gate5C_SuSiE_variant_PIP.csv.gz"))

audit <- data.table(
  metric = c(
    "LD_reference_samples", "LD_matrix_SNPs", "SuSiE_L",
    "credible_set_coverage", "WMH_signals",
    "mALPS_signals", "left_ALPS_signals", "right_ALPS_signals"
  ),
  value = c(
    489, length(snps), 10, 0.95,
    sum(!is.na(cs_tab[trait == "WMH", signal])),
    sum(!is.na(cs_tab[trait == "mALPS", signal])),
    sum(!is.na(cs_tab[trait == "left_ALPS", signal])),
    sum(!is.na(cs_tab[trait == "right_ALPS", signal]))
  )
)
fwrite(audit, file.path(fm, "Gate5C_SuSiE_audit.csv"))
writeLines(c(
  paste("R", R.version.string),
  paste("coloc", packageVersion("coloc")),
  paste("susieR", packageVersion("susieR")),
  paste("data.table", packageVersion("data.table")),
  "Random seed 20260725"
), file.path(fm, "Gate5C_SuSiE_software_versions.txt"))
print(audit)
print(coloc_tab)
