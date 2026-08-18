suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
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

project <- resolve_project_root()
g5b <- file.path(project, "Gate5B")
g5c <- file.path(project, "Gate5C")
outdir <- file.path(g5c, "01_WMH_independent")
logdir <- file.path(g5c, "07_logs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

read_gz <- function(path) fread(cmd = paste("gzip -cd", shQuote(path)))

wmh <- read_gz(file.path(outdir, "Traylor2016_WMH_chr16.tsv.gz"))
setnames(wmh, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"),
         c("SNP", "CHR", "BP", "EA", "OA", "beta", "se", "P"))
wmh[, `:=`(
  SNP = as.character(SNP),
  EA = toupper(EA),
  OA = toupper(OA),
  beta = as.numeric(beta),
  se = as.numeric(se),
  P = as.numeric(P),
  BP = as.integer(BP)
)]
wmh <- unique(wmh[is.finite(beta) & is.finite(se) & se > 0 & P > 0 & P <= 1], by = "SNP")

alps_files <- c(
  mALPS = "mALPS_Gate5B_regions.tsv.gz",
  left_ALPS = "left_ALPS_Gate5B_regions.tsv.gz",
  right_ALPS = "right_ALPS_Gate5B_regions.tsv.gz"
)

complement <- c(A = "T", T = "A", C = "G", G = "C")
summary_list <- list()
posterior_list <- list()
harmonized_list <- list()

for (phenotype in names(alps_files)) {
  a <- read_gz(file.path(g5b, "02_extracted", alps_files[[phenotype]]))
  a <- a[region_id == "ALPS_L16_011"]
  a[, `:=`(
    SNP = as.character(SNP),
    EA = toupper(EA),
    OA = toupper(OA),
    beta = as.numeric(beta),
    se = as.numeric(se),
    P = as.numeric(P),
    BP = as.integer(BP),
    EAF = as.numeric(EAF),
    N = as.numeric(N)
  )]
  a <- unique(a[is.finite(beta) & is.finite(se) & se > 0 & P > 0 & P <= 1], by = "SNP")
  m <- merge(a, wmh, by = "SNP", suffixes = c("_alps", "_wmh"))
  m[, palindromic := paste0(EA_alps, OA_alps) %chin% c("AT", "TA", "CG", "GC")]
  m[, direct := (EA_alps == EA_wmh & OA_alps == OA_wmh)]
  m[, reverse := (EA_alps == OA_wmh & OA_alps == EA_wmh)]
  m[, complement_direct := (EA_alps == complement[EA_wmh] & OA_alps == complement[OA_wmh])]
  m[, complement_reverse := (EA_alps == complement[OA_wmh] & OA_alps == complement[EA_wmh])]
  m <- m[(direct | reverse | complement_direct | complement_reverse) & !palindromic]
  m[, beta_wmh_aligned := fifelse(direct | complement_direct, beta_wmh, -beta_wmh)]
  setorder(m, SNP)
  if (nrow(m) < 50L) stop("Too few aligned SNPs for ", phenotype, ": ", nrow(m))
  if (anyDuplicated(m$SNP)) stop("Duplicate SNPs remain after harmonization for ", phenotype)

  d1 <- list(
    beta = m$beta_alps,
    varbeta = m$se_alps^2,
    snp = m$SNP,
    position = m$BP_alps,
    MAF = pmin(m$EAF, 1 - m$EAF),
    N = as.integer(round(median(m$N))),
    type = "quant"
  )
  d2 <- list(
    beta = m$beta_wmh_aligned,
    varbeta = m$se_wmh^2,
    snp = m$SNP,
    position = m$BP_wmh,
    N = 3670L,
    sdY = 1,
    type = "quant"
  )

  fit <- suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  strict <- suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-6))
  sm <- fit$summary
  vr <- as.data.table(fit$results)
  setorder(vr, -SNP.PP.H4)
  vr[, cumulative_PP_H4 := cumsum(SNP.PP.H4)]
  vr[, in_95pct_H4_credible_set := cumulative_PP_H4 <= 0.95 | shift(cumulative_PP_H4, fill = 0) < 0.95]
  vr[, `:=`(ALPS_phenotype = phenotype, external_trait = "Traylor2016_WMH_nonUKB")]

  summary_list[[phenotype]] <- data.table(
    ALPS_phenotype = phenotype,
    external_trait = "Traylor2016_WMH_nonUKB",
    sample_size_WMH = 3670L,
    common_SNPs = as.integer(sm["nsnps"]),
    PP0 = as.numeric(sm["PP.H0.abf"]),
    PP1 = as.numeric(sm["PP.H1.abf"]),
    PP2 = as.numeric(sm["PP.H2.abf"]),
    PP3 = as.numeric(sm["PP.H3.abf"]),
    PP4 = as.numeric(sm["PP.H4.abf"]),
    PP4_strict_p12_1e_6 = as.numeric(strict$summary["PP.H4.abf"]),
    top_shared_candidate_SNP = vr[1, snp],
    top_SNP_PP_H4 = vr[1, SNP.PP.H4],
    H4_95pct_credible_set_SNPs = vr[in_95pct_H4_credible_set == TRUE, .N],
    minimum_ALPS_P = min(m$P_alps),
    minimum_WMH_P = min(m$P_wmh),
    minimum_WMH_P_SNP = m[which.min(P_wmh), SNP],
    WMH_regional_P_lt_1e_6 = min(m$P_wmh) < 1e-6,
    sample_overlap = "None known; Traylor 2016 stroke cohorts do not include UK Biobank",
    model = "coloc.abf beta/varbeta; single causal variant",
    interpretation_status = ifelse(
      min(m$P_wmh) >= 1e-6,
      "LOW_POWER_NO_STRONG_WMH_REGIONAL_SIGNAL",
      ifelse(as.numeric(sm["PP.H4.abf"]) > 0.8, "INDEPENDENT_COLOC_SUPPORTED", "NO_COLOC_SUPPORT")
    )
  )
  posterior_list[[phenotype]] <- vr
  harmonized_list[[phenotype]] <- m[, .(
    ALPS_phenotype = phenotype, SNP, BP_alps, EA_alps, OA_alps, EAF,
    beta_alps, se_alps, P_alps, BP_wmh, EA_wmh, OA_wmh,
    beta_wmh, beta_wmh_aligned, se_wmh, P_wmh
  )]
}

summary_dt <- rbindlist(summary_list, fill = TRUE)
posterior_dt <- rbindlist(posterior_list, fill = TRUE)
harmonized_dt <- rbindlist(harmonized_list, fill = TRUE)
fwrite(summary_dt, file.path(outdir, "Gate5C_independent_WMH_coloc_summary.csv"))
fwrite(posterior_dt, file.path(outdir, "Gate5C_independent_WMH_coloc_variant_posteriors.csv"))
fwrite(harmonized_dt, file.path(outdir, "Gate5C_independent_WMH_harmonized.csv.gz"))

versions <- c(
  paste("R", R.version.string),
  paste("coloc", as.character(packageVersion("coloc"))),
  paste("data.table", as.character(packageVersion("data.table")))
)
writeLines(versions, file.path(logdir, "independent_WMH_coloc_software_versions.txt"))
print(summary_dt)
