#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(coloc)
  library(scales)
})

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "09_manuscript")) && dir.exists(file.path(candidate, "Gate5C"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

project <- resolve_project_root()
asset_dir <- file.path(project, "09_manuscript", "submission_assets")
fig_dir <- file.path(asset_dir, "figures")
data_dir <- file.path(asset_dir, "figure_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

harm_file <- file.path(
  project, "Gate5B", "03_coloc_results",
  "Gate5B_coloc_harmonized_regions.csv.gz"
)
susie_harm_file <- file.path(
  project, "Gate5C", "05_finemapping",
  "susie_harmonized_common.csv.gz"
)
vars_file <- file.path(
  project, "Gate5C", "05_finemapping",
  "1000G_EUR_chr16_susie.unphased.vcor1.vars"
)
ld_file <- file.path(
  project, "Gate5C", "05_finemapping",
  "1000G_EUR_chr16_susie.unphased.vcor1"
)
cs_file <- file.path(
  project, "Gate5D", "01_credible_variants",
  "chr16_SuSiE_credible_set_membership_long.csv"
)
replication_file <- file.path(
  project, "Gate5C", "01_WMH_independent",
  "Gate5C_independent_WMH_coloc_summary.csv"
)

trait_levels <- c("WMH", "Mean ALPS", "Left ALPS", "Right ALPS")
trait_map <- c(
  WMH = "WMH",
  mALPS = "Mean ALPS",
  left_ALPS = "Left ALPS",
  right_ALPS = "Right ALPS"
)
trait_colors <- c(
  "WMH" = "#A23B2A",
  "Mean ALPS" = "#1F5A94",
  "Left ALPS" = "#B56A00",
  "Right ALPS" = "#287A5B"
)
ld_colors <- c(
  "r² < 0.2" = "#4D4D4D",
  "0.2–0.4" = "#56B4E9",
  "0.4–0.6" = "#009E73",
  "0.6–0.8" = "#E69F00",
  "r² ≥ 0.8" = "#D55E00",
  "Not in LD reference" = "#C8C8C8"
)

theme_sci <- function(base_size = 8.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "#222222"),
      axis.title = element_text(color = "#222222"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(5.5, "mm"),
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle = element_text(size = rel(0.90), color = "#444444"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(fig_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(fig_dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in",
    dpi = 450, bg = "white"
  )
  tiff(
    file.path(fig_dir, paste0(stem, ".tiff")),
    width = width, height = height, units = "in", res = 600,
    compression = "lzw", type = "cairo", bg = "white"
  )
  print(plot)
  dev.off()
}

# -------------------------------------------------------------------------
# A. Stacked regional-association tracks with LD and credible-set overlays
# -------------------------------------------------------------------------

harm <- fread(cmd = paste("gzip -dc", shQuote(harm_file)))
harm <- harm[
  region_id == "ALPS_L16_011" &
    external_trait == "WMH" &
    is.finite(P_alps) & P_alps > 0 &
    is.finite(P_ext) & P_ext > 0
]

regional_alps <- harm[, .(
  SNP,
  BP = BP_alps,
  P = P_alps,
  trait = ALPS_phenotype
)]
regional_wmh <- unique(
  harm[ALPS_phenotype == "mALPS", .(
    SNP,
    BP = BP_ext,
    P = P_ext,
    trait = "WMH"
  )],
  by = "SNP"
)
regional <- rbind(regional_wmh, regional_alps, fill = TRUE)
regional[, trait_label := factor(unname(trait_map[trait]), levels = trait_levels)]
regional[, position_mb := BP / 1e6]
regional[, minus_log10_p := -log10(pmax(P, .Machine$double.xmin))]

ld_snps <- fread(vars_file, header = FALSE)[[1]]
ld <- as.matrix(fread(ld_file, header = FALSE))
storage.mode(ld) <- "double"
if (nrow(ld) != length(ld_snps) || ncol(ld) != length(ld_snps)) {
  stop("LD matrix dimensions do not match the variant list.")
}
anchor <- "rs4843555"
anchor_index <- match(anchor, ld_snps)
if (is.na(anchor_index)) stop("LD anchor rs4843555 was not found.")
ld_lookup <- data.table(
  SNP = ld_snps,
  r2_to_rs4843555 = ld[, anchor_index]^2
)
regional <- merge(regional, ld_lookup, by = "SNP", all.x = TRUE)
regional[, ld_band := cut(
  r2_to_rs4843555,
  breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
  labels = c("r² < 0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "r² ≥ 0.8"),
  right = FALSE
)]
regional[, ld_band := as.character(ld_band)]
regional[is.na(ld_band), ld_band := "Not in LD reference"]
regional[, ld_band := factor(ld_band, levels = names(ld_colors))]

cs <- fread(cs_file)
cs[, trait_label := factor(unname(trait_map[trait]), levels = trait_levels)]
regional <- merge(
  regional,
  unique(cs[, .(trait_label, SNP, in_susie_cs = TRUE)]),
  by = c("trait_label", "SNP"),
  all.x = TRUE
)
regional[is.na(in_susie_cs), in_susie_cs := FALSE]
regional[, key_snp := SNP %chin% c("rs4843550", "rs4843555")]

fwrite(
  regional[order(trait_label, BP)],
  file.path(data_dir, "Figure5A_chr16_regional_association.csv")
)

key_labels <- regional[
  key_snp == TRUE & trait == "WMH",
  .SD[which.max(minus_log10_p)],
  by = .(trait_label, SNP)
]
key_labels[, nudge_group := fifelse(SNP == "rs4843550", 1L, -1L)]

p5a <- ggplot(regional, aes(position_mb, minus_log10_p)) +
  geom_hline(
    yintercept = -log10(5e-8),
    linetype = 3, color = "#888888", linewidth = 0.30
  ) +
  geom_vline(
    xintercept = c(87231160, 87236383) / 1e6,
    linetype = 3, color = "#A0A0A0", linewidth = 0.25
  ) +
  geom_point(aes(color = ld_band), size = 0.75, alpha = 0.90) +
  geom_point(
    data = regional[in_susie_cs == TRUE],
    shape = 21, fill = NA, color = "#111111", stroke = 0.30, size = 1.45
  ) +
  geom_point(
    data = regional[key_snp == TRUE],
    shape = 23, fill = "white", color = "#111111", stroke = 0.55, size = 2.0
  ) +
  geom_text_repel(
    data = key_labels,
    aes(label = SNP),
    size = 2.05,
    min.segment.length = 0,
    segment.size = 0.25,
    nudge_x = 0.12 * key_labels$nudge_group,
    direction = "both",
    box.padding = 0.20,
    point.padding = 0.15,
    max.overlaps = Inf,
    seed = 20260728,
    show.legend = FALSE
  ) +
  facet_grid(trait_label ~ ., scales = "free_y") +
  scale_color_manual(values = ld_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(86.8, 87.7, 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Chromosome 16 position (Mb; GRCh37)",
    y = expression(-log[10](P)),
    color = expression("LD with rs4843555 (" * r^2 * ")"),
    title = "Regional association profiles in the data-selected chr16 interval",
    subtitle = "Open circles denote trait-specific SuSiE 95% credible-set variants"
  ) +
  theme_sci(7.7) +
  theme(
    legend.position = "bottom",
    panel.spacing.y = unit(1.2, "mm"),
    strip.text.y = element_text(angle = 0),
    axis.title.y = element_text(margin = margin(r = 4))
  )

# -------------------------------------------------------------------------
# B. Position-aligned SuSiE credible sets and PIPs
# -------------------------------------------------------------------------

pos_lookup <- unique(regional[, .(SNP, BP)])
cs_plot <- merge(cs, pos_lookup, by = "SNP", all.x = TRUE)
cs_plot <- cs_plot[is.finite(BP)]
cs_plot[, position_mb := BP / 1e6]
cs_plot[, trait_label := factor(unname(trait_map[trait]), levels = trait_levels)]
cs_ranges <- cs_plot[, .(
  xmin = min(position_mb),
  xmax = max(position_mb)
), by = trait_label]
lead_pip <- cs_plot[, .SD[which.max(PIP)], by = trait_label]
fwrite(
  cs_plot[order(trait_label, BP)],
  file.path(data_dir, "Figure5B_SuSiE_credible_set_positions.csv")
)

p5b <- ggplot(cs_plot, aes(position_mb, trait_label)) +
  geom_segment(
    data = cs_ranges,
    aes(x = xmin, xend = xmax, y = trait_label, yend = trait_label),
    inherit.aes = FALSE, color = "#B5B5B5", linewidth = 0.65
  ) +
  geom_point(
    aes(size = PIP, fill = trait_label),
    shape = 21, color = "#222222", stroke = 0.25, alpha = 0.95
  ) +
  geom_point(
    data = lead_pip,
    shape = 23, fill = "white", color = "#111111", stroke = 0.65, size = 2.4
  ) +
  geom_text_repel(
    data = lead_pip,
    aes(label = SNP),
    size = 2.2, nudge_x = 0.12, direction = "x",
    min.segment.length = 0, segment.size = 0.25,
    box.padding = 0.15, seed = 20260728,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = trait_colors, guide = "none") +
  scale_size_continuous(
    range = c(1.2, 5.0),
    breaks = c(0.05, 0.10, 0.20),
    name = "PIP"
  ) +
  scale_x_continuous(
    limits = range(regional$position_mb),
    breaks = seq(86.8, 87.7, 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Chromosome 16 position (Mb; GRCh37)",
    y = NULL,
    title = "SuSiE 95% credible sets"
  ) +
  theme_sci(7.7) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(face = "bold")
  )

# -------------------------------------------------------------------------
# C. Full p12 prior-sensitivity grid for classic coloc
# -------------------------------------------------------------------------

make_coloc_pair <- function(dt) {
  d1 <- list(
    beta = dt$beta_alps,
    varbeta = dt$se_alps^2,
    snp = dt$SNP,
    position = dt$BP_alps,
    MAF = pmin(dt$EAF_alps, 1 - dt$EAF_alps),
    N = 31021L,
    type = "quant"
  )
  d2 <- list(
    beta = dt$beta_ext,
    varbeta = dt$se_ext^2,
    snp = dt$SNP,
    position = dt$BP_ext,
    MAF = pmin(dt$EAF_ext, 1 - dt$EAF_ext),
    N = 48454L,
    type = "quant"
  )
  list(d1 = d1, d2 = d2)
}

p12_grid <- c(1e-7, 1e-6, 1e-5, 1e-4)
sensitivity <- rbindlist(lapply(
  c("mALPS", "left_ALPS", "right_ALPS"),
  function(pheno) {
    dt <- harm[
      ALPS_phenotype == pheno &
        is.finite(beta_alps) & is.finite(se_alps) & se_alps > 0 &
        is.finite(beta_ext) & is.finite(se_ext) & se_ext > 0 &
        is.finite(EAF_alps) & EAF_alps > 0 & EAF_alps < 1 &
        is.finite(EAF_ext) & EAF_ext > 0 & EAF_ext < 1
    ]
    pair <- make_coloc_pair(dt)
    rbindlist(lapply(p12_grid, function(p12_value) {
      fit <- suppressWarnings(coloc.abf(
        dataset1 = pair$d1,
        dataset2 = pair$d2,
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = p12_value
      ))
      data.table(
        ALPS_phenotype = pheno,
        p12 = p12_value,
        PP0 = as.numeric(fit$summary["PP.H0.abf"]),
        PP1 = as.numeric(fit$summary["PP.H1.abf"]),
        PP2 = as.numeric(fit$summary["PP.H2.abf"]),
        PP3 = as.numeric(fit$summary["PP.H3.abf"]),
        PP4 = as.numeric(fit$summary["PP.H4.abf"]),
        common_SNPs = as.integer(fit$summary["nsnps"])
      )
    }))
  }
))
sensitivity[, trait_label := factor(
  unname(trait_map[ALPS_phenotype]),
  levels = trait_levels[-1]
)]
sensitivity_long <- melt(
  sensitivity,
  id.vars = c("ALPS_phenotype", "trait_label", "p12", "common_SNPs"),
  measure.vars = c("PP3", "PP4"),
  variable.name = "hypothesis",
  value.name = "posterior"
)
fwrite(
  sensitivity_long[order(trait_label, p12, hypothesis)],
  file.path(data_dir, "Figure5C_coloc_prior_sensitivity_full.csv")
)

p5c <- ggplot(
  sensitivity_long,
  aes(p12, posterior, color = trait_label, linetype = hypothesis, group = hypothesis)
) +
  geom_hline(yintercept = 0.80, linetype = 3, color = "#777777", linewidth = 0.35) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.35) +
  facet_wrap(~ trait_label, nrow = 1) +
  scale_x_log10(
    breaks = p12_grid,
    labels = c(expression(10^-7), expression(10^-6), expression(10^-5), expression(10^-4))
  ) +
  scale_color_manual(values = trait_colors[trait_levels[-1]], guide = "none") +
  scale_linetype_manual(values = c(PP3 = 2, PP4 = 1)) +
  scale_y_continuous(limits = c(0, 1.01), breaks = seq(0, 1, 0.2)) +
  labs(
    x = expression("Shared-association prior " * p[12]),
    y = "Posterior probability",
    linetype = NULL,
    title = expression("Conventional-coloc sensitivity across " * p[12] * " priors")
  ) +
  theme_sci(7.6) +
  theme(
    legend.position = "bottom",
    panel.spacing.x = unit(2.2, "mm")
  )

# -------------------------------------------------------------------------
# D. Primary and non-UKB formal evaluation
# -------------------------------------------------------------------------

primary <- sensitivity[p12 == 1e-5, .(
  ALPS_phenotype,
  trait_label,
  dataset = "Primary WMH (overlapping)",
  PP3,
  PP4
)]
replication <- fread(replication_file)
replication[, trait_label := factor(
  unname(trait_map[ALPS_phenotype]),
  levels = trait_levels[-1]
)]
external <- replication[, .(
  ALPS_phenotype,
  trait_label,
  dataset = "Traylor WMH (non-UKB)",
  PP3,
  PP4
)]
evaluation <- rbind(primary, external)
evaluation_long <- melt(
  evaluation,
  id.vars = c("ALPS_phenotype", "trait_label", "dataset"),
  measure.vars = c("PP3", "PP4"),
  variable.name = "hypothesis",
  value.name = "posterior"
)
fwrite(
  evaluation_long[order(dataset, trait_label, hypothesis)],
  file.path(data_dir, "Figure5D_primary_external_evaluation.csv")
)

p5d <- ggplot(
  evaluation_long,
  aes(trait_label, posterior, fill = hypothesis)
) +
  geom_hline(yintercept = 0.80, linetype = 3, color = "#777777", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  facet_wrap(~ dataset, nrow = 1, labeller = label_wrap_gen(width = 18)) +
  scale_fill_manual(values = c(PP3 = "#BDBDBD", PP4 = "#3D6594")) +
  scale_y_continuous(limits = c(0, 1.01), breaks = seq(0, 1, 0.2)) +
  labs(
    x = NULL,
    y = "Posterior probability",
    fill = NULL,
    title = "Formal regional evaluation by WMH dataset",
    subtitle = "CHARGE: P-value-truncated lookup; PP4 not estimable"
  ) +
  theme_sci(7.6) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

figure5 <- p5a / p5b / (p5c | p5d) +
  plot_layout(heights = c(3.9, 1.65, 2.2), widths = c(1.15, 0.85)) +
  plot_annotation(
    title = "Data-selected chr16 ALPS–WMH regional evaluation",
    subtitle = paste(
      "The region was selected after the primary WMH screen;",
      "regional analyses are internal and susceptible to winner selection."
    ),
    tag_levels = "A",
    caption = paste0(
      "Primary ALPS and WMH GWAS substantially overlap. SuSiE addresses multiple ",
      "association signals but not cross-trait sampling covariance from overlap.\n",
      "rs4843550 was retained in regional plots but excluded from the LD-based ",
      "SuSiE set because its palindromic allele orientation was unresolved."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 10.2, hjust = 0),
      plot.subtitle = element_text(size = 7.8, hjust = 0),
      plot.tag = element_text(face = "bold", size = 10),
      plot.caption = element_text(
        size = 7.0, hjust = 0, color = "#222222", lineheight = 1.05
      )
    )
  )

save_figure(
  figure5,
  "Figure5_chr16_colocalization_finemapping_replication",
  width = 7.2,
  height = 9.2
)

message("Figure 5 Q2 redesign completed.")
