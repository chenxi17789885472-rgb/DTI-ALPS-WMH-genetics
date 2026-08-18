suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(jsonlite)
})

set.seed(20260726)

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5C", "10_mALPS_to_WMH_MR"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

project <- resolve_project_root()
stage <- file.path(project, "Gate5C", "10_mALPS_to_WMH_MR")
harm_file <- file.path(stage, "02_harmonization", "mALPS_Traylor2016_WMH_harmonized.tsv")
result_dir <- file.path(stage, "03_results")
sens_dir <- file.path(stage, "04_sensitivity")
log_dir <- file.path(stage, "05_logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

x <- fread(harm_file)
required_columns <- c(
  "SNP", "beta_exposure", "se_exposure", "EA_exposure", "OA_exposure",
  "EAF_exposure", "p_exposure", "N_exposure", "beta_outcome_aligned",
  "se_outcome", "p_outcome", "N_outcome", "F_exposure", "mr_keep", "proxy_used"
)
missing_columns <- setdiff(required_columns, names(x))
if (length(missing_columns)) stop("Harmonized MR input is missing: ", paste(missing_columns, collapse = ", "))
if (nrow(x) < 10L) stop("Fewer than 10 harmonized instruments: ", nrow(x))
if (anyDuplicated(x$SNP)) stop("Duplicate instruments in harmonized MR input.")
if (any(!is.finite(x$beta_exposure)) || any(!is.finite(x$se_exposure)) ||
    any(!is.finite(x$beta_outcome_aligned)) || any(!is.finite(x$se_outcome))) {
  stop("Non-finite beta or standard-error values in harmonized MR input.")
}
if (any(x$se_exposure <= 0) || any(x$se_outcome <= 0)) stop("Non-positive standard errors in MR input.")
if (any(x$F_exposure <= 10) || !all(x$mr_keep %in% TRUE)) stop("Instrument-strength or mr_keep gate failed.")

dat <- data.table(
  SNP = x$SNP,
  beta.exposure = x$beta_exposure,
  se.exposure = x$se_exposure,
  effect_allele.exposure = x$EA_exposure,
  other_allele.exposure = x$OA_exposure,
  eaf.exposure = x$EAF_exposure,
  pval.exposure = x$p_exposure,
  samplesize.exposure = x$N_exposure,
  beta.outcome = x$beta_outcome_aligned,
  se.outcome = x$se_outcome,
  effect_allele.outcome = x$EA_exposure,
  other_allele.outcome = x$OA_exposure,
  eaf.outcome = x$EAF_exposure,
  pval.outcome = x$p_outcome,
  samplesize.outcome = x$N_outcome,
  exposure = "Mean DTI-ALPS index",
  outcome = "WMH burden (Traylor 2016 non-UKB)",
  id.exposure = "mALPS_Figshare_25515331",
  id.outcome = "Traylor2016_WMH_nonUKB",
  mr_keep = TRUE,
  stringsAsFactors = FALSE
)

parameters <- default_parameters()
parameters$nboot <- 10000
mr_res <- as.data.table(mr(
  dat,
  parameters = parameters,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_raps")
))
mr_res[, `:=`(
  ci_lower = b - qnorm(0.975) * se,
  ci_upper = b + qnorm(0.975) * se,
  direction = fifelse(b < 0, "negative", fifelse(b > 0, "positive", "null")),
  analysis_status = "exploratory_nonUKB_low_power"
)]
fwrite(mr_res, file.path(result_dir, "mALPS_to_Traylor2016_WMH_MR_results.csv"))

het <- as.data.table(mr_heterogeneity(
  dat, method_list = c("mr_ivw", "mr_egger_regression")
))
fwrite(het, file.path(sens_dir, "Cochran_Q_results.csv"))

pleio <- as.data.table(mr_pleiotropy_test(dat))
fwrite(pleio, file.path(sens_dir, "MR_Egger_intercept.csv"))

loo <- as.data.table(mr_leaveoneout(dat, method = mr_ivw))
loo[, `:=`(
  ci_lower = b - qnorm(0.975) * se,
  ci_upper = b + qnorm(0.975) * se,
  direction = fifelse(b < 0, "negative", fifelse(b > 0, "positive", "null"))
)]
fwrite(loo, file.path(sens_dir, "leave_one_out_results.csv"))

single <- as.data.table(mr_singlesnp(dat, all_method = c("mr_ivw")))
fwrite(single, file.path(sens_dir, "single_SNP_results.csv"))

run_subset <- function(d, label) {
  z <- as.data.table(mr(
    d,
    parameters = parameters,
    method_list = c("mr_ivw", "mr_weighted_median")
  ))
  z[, `:=`(
    subset = label,
    ci_lower = b - qnorm(0.975) * se,
    ci_upper = b + qnorm(0.975) * se
  )]
  z[]
}

proxy_snps <- x[proxy_used == TRUE, SNP]
subset_results <- rbindlist(list(
  run_subset(dat, "all_10_IV"),
  run_subset(dat[!SNP %in% proxy_snps], "direct_only_8_IV"),
  run_subset(dat[SNP != "rs4843555"], "chr16_IV_excluded"),
  run_subset(dat[SNP != "rs10817103"], "strongest_WMH_SNP_excluded")
), fill = TRUE)
fwrite(subset_results, file.path(sens_dir, "prespecified_subset_sensitivity.csv"))

presso_data <- as.data.frame(dat)
rownames(presso_data) <- presso_data$SNP
presso <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  data = presso_data,
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  SignifThreshold = 0.05,
  NbDistribution = 10000,
  seed = 20260726
)
saveRDS(presso, file.path(sens_dir, "MR_PRESSO_full_result.rds"))

global <- presso[["MR-PRESSO results"]][["Global Test"]]
outlier <- presso[["MR-PRESSO results"]][["Outlier Test"]]
distortion <- presso[["MR-PRESSO results"]][["Distortion Test"]]
main_presso <- as.data.table(presso[["Main MR results"]])
fwrite(main_presso, file.path(sens_dir, "MR_PRESSO_main_MR.csv"))
if (!is.null(outlier)) {
  fwrite(as.data.table(outlier, keep.rownames = "SNP"),
         file.path(sens_dir, "MR_PRESSO_outlier_test.csv"))
}

parse_p <- function(v) suppressWarnings(as.numeric(sub("^<", "", as.character(v))))
outlier_snps <- character()
if (!is.null(outlier)) {
  odt <- as.data.table(outlier, keep.rownames = "SNP")
  outlier_snps <- odt[parse_p(Pvalue) < 0.05, SNP]
}

presso_summary <- data.table(
  nsnp = nrow(dat),
  global_RSSobs = as.numeric(global$RSSobs),
  global_P_display = as.character(global$Pvalue),
  global_P_numeric_bound = parse_p(global$Pvalue),
  outlier_count = length(outlier_snps),
  outlier_SNPs = paste(outlier_snps, collapse = ";"),
  distortion_coefficient = if (is.null(distortion)) NA_real_ else
    as.numeric(distortion[["Distortion Coefficient"]]),
  distortion_P = if (is.null(distortion)) NA_real_ else
    parse_p(distortion$Pvalue)
)
fwrite(presso_summary, file.path(sens_dir, "MR_PRESSO_summary.csv"))

# Approximate Steiger directionality for standardized quantitative phenotypes.
steiger <- x[, .(
  SNP,
  r2_exposure = 2 * EAF_exposure * (1 - EAF_exposure) * beta_exposure^2,
  r2_outcome = 2 * EAF_exposure * (1 - EAF_exposure) * beta_outcome_aligned^2
)]
steiger[, correct_per_SNP := r2_exposure > r2_outcome]
fwrite(steiger, file.path(sens_dir, "Steiger_per_SNP_approx.csv"))
steiger_summary <- steiger[, .(
  nsnp = .N,
  sum_r2_exposure = sum(r2_exposure),
  sum_r2_outcome = sum(r2_outcome),
  correct_direction = sum(r2_exposure) > sum(r2_outcome),
  SNP_fraction_correct = mean(correct_per_SNP),
  caveat = "Approximate: both traits treated as standardized quantitative phenotypes; outcome EAF proxied by exposure EAF."
)]
fwrite(steiger_summary, file.path(sens_dir, "Steiger_summary_approx.csv"))

audit <- list(
  exposure_IV_original = 12L,
  direct_outcome_matches = 10L,
  conservative_nonpalindromic_direct = 8L,
  proxies_added = sum(x$proxy_used),
  final_IV = nrow(x),
  minimum_F = min(x$F_exposure),
  median_F = median(x$F_exposure),
  go_threshold = 10L,
  gate_status = ifelse(nrow(x) >= 10L, "GO_EXPLORATORY_MR", "NO_GO"),
  exposure_sample_median = median(x$N_exposure),
  outcome_sample_size = 3670L,
  sample_overlap = "None known",
  analysis_limit = "Small stroke-cohort WMH outcome; low power; exactly at IV threshold."
)
write_json(audit, file.path(log_dir, "mALPS_WMH_MR_audit.json"),
           pretty = TRUE, auto_unbox = TRUE)

writeLines(c(
  paste("R", R.version.string),
  paste("TwoSampleMR", packageVersion("TwoSampleMR")),
  paste("MRPRESSO", packageVersion("MRPRESSO")),
  paste("data.table", packageVersion("data.table"))
), file.path(log_dir, "software_versions.txt"))

print(mr_res)
print(het)
print(pleio)
print(presso_summary)
print(steiger_summary)
