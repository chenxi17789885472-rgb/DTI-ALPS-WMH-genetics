suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(ggplot2)
})

set.seed(20260725)

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5B"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

project <- resolve_project_root()
gate <- file.path(project, "Gate5B")
input_dir <- file.path(gate, "02_extracted")
result_dir <- file.path(gate, "03_coloc_results")
figure_dir <- file.path(gate, "04_figures")
log_dir <- file.path(gate, "05_logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

read_gz <- function(path) fread(cmd = paste("gzip -cd", shQuote(path)))

require_columns <- function(dat, columns, label) {
  missing <- setdiff(columns, names(dat))
  if (length(missing)) stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
}

alps <- list(
  mALPS = read_gz(file.path(input_dir, "mALPS_Gate5B_regions.tsv.gz")),
  left_ALPS = read_gz(file.path(input_dir, "left_ALPS_Gate5B_regions.tsv.gz")),
  right_ALPS = read_gz(file.path(input_dir, "right_ALPS_Gate5B_regions.tsv.gz"))
)
external <- list(
  WMH = read_gz(file.path(input_dir, "WMH_chr16.tsv.gz")),
  sleep_duration = read_gz(file.path(input_dir, "sleep_duration_chr17.tsv.gz")),
  insomnia = read_gz(file.path(input_dir, "insomnia_chr17.tsv.gz"))
)

specs <- data.table(
  external_trait = c("WMH", "sleep_duration", "insomnia"),
  region_id = c("ALPS_L16_011", "ALPS_L17_012", "ALPS_L17_012"),
  type = c("quant", "quant", "cc"),
  N = c(48454L, 446118L, 237627L),
  case_fraction = c(NA_real_, NA_real_, 129270 / 237627),
  sample_overlap = c(
    "Possible major UKB imaging overlap (WMH includes 26,788 UKB participants)",
    "UKB overlap likely; up to approximately 31,021 ALPS participants",
    "UKB overlap likely; up to approximately 31,021 ALPS participants"
  )
)

complement <- c(A = "T", T = "A", C = "G", G = "C")

clean_dataset <- function(dat, apply_info = FALSE) {
  dat <- copy(dat)
  require_columns(
    dat,
    c("SNP", "EA", "OA", "EAF", "INFO", "beta", "se", "P", "BP", "N"),
    "Regional summary-statistics input"
  )
  dat[, `:=`(
    SNP = as.character(SNP),
    EA = toupper(as.character(EA)),
    OA = toupper(as.character(OA)),
    EAF = as.numeric(EAF),
    INFO = as.numeric(INFO),
    beta = as.numeric(beta),
    se = as.numeric(se),
    P = as.numeric(P),
    BP = as.integer(BP),
    N = as.numeric(N)
  )]
  dat <- dat[
    SNP != "" &
      EA %chin% names(complement) &
      OA %chin% names(complement) &
      is.finite(EAF) & EAF > 0 & EAF < 1 &
      is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(P) & P > 0 & P <= 1
  ]
  if (apply_info) {
    dat <- dat[is.finite(INFO) & INFO >= 0.8]
  }
  setorder(dat, P)
  unique(dat, by = "SNP")
}

harmonize_pair <- function(a, e) {
  a <- clean_dataset(a, apply_info = FALSE)
  e <- clean_dataset(e, apply_info = any(is.finite(suppressWarnings(as.numeric(e$INFO)))))
  m <- merge(a, e, by = "SNP", suffixes = c("_alps", "_ext"))
  m[, direct_set := (
    (EA_alps == EA_ext & OA_alps == OA_ext) |
      (EA_alps == OA_ext & OA_alps == EA_ext)
  )]
  m[, complement_set := (
    (EA_alps == complement[EA_ext] & OA_alps == complement[OA_ext]) |
      (EA_alps == complement[OA_ext] & OA_alps == complement[EA_ext])
  )]
  m[, palindromic := paste0(EA_alps, OA_alps) %chin% c("AT", "TA", "CG", "GC")]
  m[, allele_status := fifelse(
    direct_set & palindromic,
    "palindromic_rsid_match",
    fifelse(
      direct_set,
      "direct_or_reversed",
      fifelse(complement_set, "strand_complement", "allele_mismatch")
    )
  )]
  m[, external_alignment_sign := fifelse(
    palindromic,
    NA_real_,
    fifelse(
      EA_ext == EA_alps | complement[EA_ext] == EA_alps,
      1,
      fifelse(OA_ext == EA_alps | complement[OA_ext] == EA_alps, -1, NA_real_)
    )
  )]
  m <- m[direct_set | complement_set]
  m[, beta_ext_aligned := beta_ext * external_alignment_sign]
  setorder(m, SNP)
  m
}

results <- list()
variant_results <- list()
allele_qc <- list()
regional_qc <- list()

for (i in seq_len(nrow(specs))) {
  spec <- specs[i]
  region <- spec$region_id
  ext_key <- spec$external_trait
  for (alps_key in names(alps)) {
    pair_id <- paste(region, alps_key, ext_key, sep = "__")
    a0 <- alps[[alps_key]][region_id == region]
    e0 <- external[[ext_key]][region_id == region]
    m <- harmonize_pair(a0, e0)

    if (!nrow(m)) {
      results[[pair_id]] <- data.table(
        pair_id = pair_id, region_id = region, ALPS_phenotype = alps_key,
        external_trait = ext_key, common_SNPs = 0L, status = "NO_ALLELE_ALIGNED_SNPS"
      )
      next
    }

    allele_qc[[pair_id]] <- data.table(
      pair_id = pair_id,
      region_id = region,
      ALPS_phenotype = alps_key,
      external_trait = ext_key,
      raw_ALPS_rows = nrow(a0),
      raw_external_rows = nrow(e0),
      common_rsid_before_allele_QC = length(intersect(a0$SNP, e0$SNP)),
      common_SNPs_after_allele_QC = nrow(m),
      palindromic_SNPs = sum(m$palindromic),
      strand_complement_SNPs = sum(m$allele_status == "strand_complement"),
      minimum_ALPS_P = min(m$P_alps),
      minimum_external_P = min(m$P_ext)
    )

    if (nrow(m) < 50) {
      results[[pair_id]] <- data.table(
        pair_id = pair_id,
        region_id = region,
        ALPS_phenotype = alps_key,
        external_trait = ext_key,
        common_SNPs = nrow(m),
        status = "INSUFFICIENT_COMMON_SNPS"
      )
      next
    }

    d_alps <- list(
      beta = m$beta_alps,
      varbeta = m$se_alps^2,
      snp = m$SNP,
      position = m$BP_alps,
      MAF = pmin(m$EAF_alps, 1 - m$EAF_alps),
      N = as.integer(round(median(m$N_alps))),
      type = "quant"
    )

    if (ext_key == "insomnia") {
      d_ext <- list(
        pvalues = m$P_ext,
        snp = m$SNP,
        position = m$BP_ext,
        MAF = pmin(m$EAF_ext, 1 - m$EAF_ext),
        N = spec$N,
        type = "cc",
        s = spec$case_fraction
      )
    } else {
      d_ext <- list(
        beta = m$beta_ext,
        varbeta = m$se_ext^2,
        snp = m$SNP,
        position = m$BP_ext,
        MAF = pmin(m$EAF_ext, 1 - m$EAF_ext),
        N = spec$N,
        type = "quant"
      )
    }

    check_dataset(d_alps, warn.minp = 1e-4)
    check_dataset(d_ext, warn.minp = 1e-4)
    fit <- suppressWarnings(coloc.abf(
      dataset1 = d_alps,
      dataset2 = d_ext,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = 1e-5
    ))
    strict_fit <- suppressWarnings(coloc.abf(
      dataset1 = d_alps,
      dataset2 = d_ext,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = 1e-6
    ))

    insomnia_logor_PP4 <- NA_real_
    if (ext_key == "insomnia") {
      scale_factor <- spec$case_fraction * (1 - spec$case_fraction)
      d_ext_logor <- list(
        beta = m$beta_ext / scale_factor,
        varbeta = (m$se_ext / scale_factor)^2,
        snp = m$SNP,
        position = m$BP_ext,
        MAF = pmin(m$EAF_ext, 1 - m$EAF_ext),
        N = spec$N,
        type = "cc",
        s = spec$case_fraction
      )
      logor_fit <- suppressWarnings(coloc.abf(
        dataset1 = d_alps,
        dataset2 = d_ext_logor,
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = 1e-5
      ))
      insomnia_logor_PP4 <- as.numeric(logor_fit$summary["PP.H4.abf"])
    }

    sm <- fit$summary
    ss <- strict_fit$summary
    vr <- as.data.table(fit$results)
    setorder(vr, -SNP.PP.H4)
    vr[, cumulative_PP_H4 := cumsum(SNP.PP.H4)]
    vr[, in_95pct_H4_credible_set := cumulative_PP_H4 <= 0.95 | shift(cumulative_PP_H4, fill = 0) < 0.95]
    top <- vr[1]
    credible <- vr[in_95pct_H4_credible_set == TRUE]

    results[[pair_id]] <- data.table(
      pair_id = pair_id,
      region_id = region,
      ALPS_phenotype = alps_key,
      external_trait = ext_key,
      common_SNPs = as.integer(sm["nsnps"]),
      PP0 = as.numeric(sm["PP.H0.abf"]),
      PP1 = as.numeric(sm["PP.H1.abf"]),
      PP2 = as.numeric(sm["PP.H2.abf"]),
      PP3 = as.numeric(sm["PP.H3.abf"]),
      PP4 = as.numeric(sm["PP.H4.abf"]),
      PP4_strict_p12_1e_6 = as.numeric(ss["PP.H4.abf"]),
      insomnia_logOR_approx_PP4 = insomnia_logor_PP4,
      PP4_gt_0_8 = as.numeric(sm["PP.H4.abf"]) > 0.8,
      top_shared_candidate_SNP = top$snp,
      top_SNP_PP_H4 = top$SNP.PP.H4,
      H4_95pct_credible_set_SNPs = nrow(credible),
      minimum_ALPS_P = min(m$P_alps),
      minimum_external_P = min(m$P_ext),
      both_traits_regional_P_lt_1e_6 = min(m$P_alps) < 1e-6 & min(m$P_ext) < 1e-6,
      sample_overlap = spec$sample_overlap,
      overlap_sensitive = TRUE,
      status = "OK",
      prior_p1 = 1e-4,
      prior_p2 = 1e-4,
      prior_p12 = 1e-5,
      model = ifelse(ext_key == "insomnia", "coloc.abf; insomnia p-value mode", "coloc.abf; beta/varbeta mode"),
      assumption = "At most one causal variant per trait in region"
    )
    vr[, `:=`(
      pair_id = pair_id,
      region_id = region,
      ALPS_phenotype = alps_key,
      external_trait = ext_key
    )]
    variant_results[[pair_id]] <- vr

    regional_qc[[pair_id]] <- m[, .(
      pair_id = pair_id,
      region_id = region,
      ALPS_phenotype = alps_key,
      external_trait = ext_key,
      SNP,
      BP_alps,
      EA_alps,
      OA_alps,
      EAF_alps,
      beta_alps,
      se_alps,
      P_alps,
      BP_ext,
      EA_ext,
      OA_ext,
      EAF_ext,
      beta_ext,
      beta_ext_aligned,
      se_ext,
      P_ext,
      INFO_ext,
      allele_status,
      palindromic
    )]
  }
}

summary_dt <- rbindlist(results, fill = TRUE)
variant_dt <- rbindlist(variant_results, fill = TRUE)
allele_qc_dt <- rbindlist(allele_qc, fill = TRUE)
regional_qc_dt <- rbindlist(regional_qc, fill = TRUE)
setorder(summary_dt, -PP4)

fwrite(summary_dt, file.path(result_dir, "Gate5B_coloc_summary.csv"))
fwrite(variant_dt, file.path(result_dir, "Gate5B_coloc_variant_posteriors.csv"))
fwrite(allele_qc_dt, file.path(result_dir, "Gate5B_coloc_allele_QC.csv"))
fwrite(regional_qc_dt, file.path(result_dir, "Gate5B_coloc_harmonized_regions.csv.gz"))

gate_summary <- summary_dt[, .(
  tested_pairs = .N,
  PP4_gt_0_8_pairs = sum(PP4_gt_0_8, na.rm = TRUE),
  maximum_PP4 = max(PP4, na.rm = TRUE),
  maximum_strict_PP4 = max(PP4_strict_p12_1e_6, na.rm = TRUE),
  all_pairs_overlap_sensitive = all(overlap_sensitive)
)]
fwrite(gate_summary, file.path(result_dir, "Gate5B_coloc_gate_summary.csv"))

p <- ggplot(
  summary_dt[status == "OK"],
  aes(x = interaction(external_trait, ALPS_phenotype, sep = " × "), y = PP4, fill = external_trait)
) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#64748B") +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c(WMH = "#7A5195", sleep_duration = "#2F4B7C", insomnia = "#EF8354")) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", color = "#17365D"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = "Gate 5B regional colocalization",
    subtitle = "Dashed line: PP4 = 0.80; all estimates are sample-overlap sensitive",
    x = NULL,
    y = "Posterior probability of a shared variant (PP4)",
    fill = NULL
  )
ggsave(file.path(figure_dir, "Gate5B_coloc_PP4.png"), p, width = 7.4, height = 5.4, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "Gate5B_coloc_PP4.pdf"), p, width = 7.4, height = 5.4, device = cairo_pdf, bg = "white")

writeLines(
  c(
    paste("R version:", R.version.string),
    paste("coloc version:", as.character(packageVersion("coloc"))),
    paste("data.table version:", as.character(packageVersion("data.table"))),
    paste("ggplot2 version:", as.character(packageVersion("ggplot2"))),
    paste("Completed:", Sys.time())
  ),
  file.path(log_dir, "software_versions.txt")
)

print(summary_dt)
