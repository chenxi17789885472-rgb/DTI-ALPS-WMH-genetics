#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
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
table_dir <- file.path(asset_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c(
  blue = "#0072B2", orange = "#E69F00", green = "#009E73",
  vermilion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
  yellow = "#F0E442", black = "#222222", grey = "#777777",
  lightgrey = "#D9D9D9", verylight = "#F3F5F7"
)

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(color = cols["black"]),
      axis.text = element_text(color = cols["black"]),
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      strip.background = element_rect(fill = "white", color = cols["lightgrey"]),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(3.5, "mm"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

save_figure <- function(plot, stem, width, height) {
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  tif_path <- file.path(fig_dir, paste0(stem, ".tiff"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  ggsave(png_path, plot, width = width, height = height, units = "in",
         dpi = 300, bg = "white")
  tiff(
    tif_path, width = width, height = height, units = "in", res = 600,
    compression = "lzw", type = "cairo", bg = "white"
  )
  print(plot)
  dev.off()
}

pretty_trait <- function(x) {
  x <- gsub("_", " ", x)
  x <- sub("^mALPS$", "Mean ALPS", x)
  x <- sub("^left ALPS$", "Left ALPS", x)
  x <- sub("^right ALPS$", "Right ALPS", x)
  x <- sub("^left ALPS$", "Left ALPS", x)
  x <- sub("^right ALPS$", "Right ALPS", x)
  x
}

# ---------------------------------------------------------------------------
# Figure 1: study design schematic
# ---------------------------------------------------------------------------

figure1_svg <- file.path(
  fig_dir, "diagram", "Figure1_study_design_redesign.svg"
)
if (!file.exists(figure1_svg)) {
nodes <- data.table(
  id = 1:8,
  x = c(1, 3, 5, 7, 7, 5, 3, 1),
  y = c(3, 3, 3, 3, 1, 1, 1, 1),
  label = c(
    "UK Biobank\nDTI-ALPS GWAS\nmean · left · right",
    "Genetic architecture\nloci · LD clumping\nSNP heritability",
    "Targeted cross-trait screen\n13 prespecified\nbrain-health traits",
    "Genome-wide sharing\nLDSC genetic\ncorrelation",
    "Local sharing\nchr16 coloc\nand SuSiE",
    "Non-UKB evidence\nTraylor coloc\nCHARGE screen",
    "Functional annotation\nVEP · GTEx · MAGMA\npathways",
    "Exploratory MR\nmALPS to WMH\nnon-UKB outcome"
  ),
  fill = c(cols["blue"], cols["sky"], cols["orange"], cols["green"],
           cols["purple"], cols["vermilion"], "#8DA0CB", "#66C2A5")
)

edges <- data.table(
  x = c(1.55, 3.55, 5.55, 7, 6.45, 4.45, 2.45),
  y = c(3, 3, 3, 2.65, 1, 1, 1),
  xend = c(2.45, 4.45, 6.45, 7, 5.55, 3.55, 1.55),
  yend = c(3, 3, 3, 1.35, 1, 1, 1)
)

fig1 <- ggplot() +
  geom_rect(
    data = nodes,
    aes(xmin = x - 0.72, xmax = x + 0.72, ymin = y - 0.42, ymax = y + 0.42,
        fill = I(fill)),
    color = "white", linewidth = 0.7
  ) +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.7, color = cols["grey"],
    arrow = arrow(length = unit(2.5, "mm"), type = "closed")
  ) +
  geom_text(
    data = nodes, aes(x = x, y = y, label = label),
    color = "white", size = 2.55, lineheight = 0.95, fontface = "bold"
  ) +
  annotate(
    "text", x = 4, y = 3.72,
    label = "Discovery and prioritization", fontface = "bold", size = 3.6,
    color = cols["black"]
  ) +
  annotate(
    "text", x = 4, y = 0.28,
    label = "Validation, annotation, and directional analysis",
    fontface = "bold", size = 3.6, color = cols["black"]
  ) +
  coord_cartesian(xlim = c(0.25, 7.75), ylim = c(0.05, 3.95), clip = "off") +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(8, 8, 8, 8))

fwrite(nodes[, .(step = id, label)], file.path(data_dir, "Figure1_workflow_nodes.csv"))
save_figure(fig1, "Figure1_study_design", 7.2, 4.4)
} else {
  message(
    "Figure 1 is maintained as a publication SVG: ",
    figure1_svg,
    ". The R fallback was not used."
  )
}

# ---------------------------------------------------------------------------
# Figure 2: ALPS genetic architecture
# ---------------------------------------------------------------------------

qc <- fread(file.path(project, "01_ALPS_GWAS", "QC", "ALPS_input_QC.tsv"))
qc[, phenotype_label := factor(
  fifelse(phenotype == "mALPS", "Mean ALPS",
          fifelse(phenotype == "left_ALPS", "Left ALPS", "Right ALPS")),
  levels = c("Mean ALPS", "Left ALPS", "Right ALPS")
)]
count_long <- melt(
  qc[, .(phenotype_label, gws_raw_rows, gws_biallelic_unique_candidates,
         strict_independent_sentinels)],
  id.vars = "phenotype_label", variable.name = "stage", value.name = "count"
)
count_long[, stage := factor(stage,
  levels = c("gws_raw_rows", "gws_biallelic_unique_candidates",
             "strict_independent_sentinels"),
  labels = c("Genome-wide\nsignificant", "QC-retained", "Independent\nsentinels")
)]
fwrite(count_long, file.path(data_dir, "Figure2A_variant_counts.csv"))

p2a <- ggplot(count_long, aes(phenotype_label, count, fill = stage)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  geom_text(
    aes(label = count), position = position_dodge(width = 0.76),
    vjust = -0.3, size = 2.6
  ) +
  scale_fill_manual(values = unname(c(cols["lightgrey"], cols["sky"], cols["blue"]))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Variant count", fill = NULL, title = "Variant prioritization") +
  theme_pub() +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold"))

h2 <- fread(file.path(project, "03_LDSC", "h2", "ALPS_LDSC_heritability_results.csv"))
h2[, phenotype_label := factor(
  fifelse(phenotype == "mALPS", "Mean ALPS",
          fifelse(phenotype == "left_ALPS", "Left ALPS", "Right ALPS")),
  levels = rev(c("Mean ALPS", "Left ALPS", "Right ALPS"))
)]
fwrite(h2[, .(phenotype, h2_observed, SE, CI_lower, CI_upper, P, LDSC_intercept,
              attenuation_ratio, regression_SNPs)],
       file.path(data_dir, "Figure2B_heritability.csv"))

p2b <- ggplot(h2, aes(h2_observed, phenotype_label)) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.12,
                orientation = "y", linewidth = 0.7, color = cols["blue"]) +
  geom_point(size = 3.0, color = cols["blue"]) +
  geom_text(aes(label = sprintf("%.3f", h2_observed)), nudge_y = 0.25,
            size = 2.7, color = cols["black"]) +
  scale_x_continuous(limits = c(0.20, 0.37), breaks = seq(0.20, 0.35, 0.05)) +
  labs(x = expression("SNP heritability (" * h[SNP]^2 * ")"),
       y = NULL, title = "Common-variant heritability") +
  theme_pub()

loci <- fread(file.path(project, "02_genetic_architecture", "loci",
                        "ALPS_merged_loci.tsv"))
membership <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  phenos <- strsplit(loci$phenotypes[i], ";", fixed = TRUE)[[1]]
  data.table(
    locus_id = loci$locus_id[i],
    lead_SNP = loci$lead_SNP[i],
    CHR = loci$CHR[i],
    phenotype = c("mALPS", "left_ALPS", "right_ALPS"),
    present = as.integer(c("mALPS", "left_ALPS", "right_ALPS") %in% phenos),
    phenotype_count = loci$phenotype_count[i]
  )
}))
membership[, locus_label := paste0("chr", CHR, "\n", lead_SNP)]
membership[, phenotype_label := factor(
  fifelse(phenotype == "mALPS", "Mean ALPS",
          fifelse(phenotype == "left_ALPS", "Left ALPS", "Right ALPS")),
  levels = rev(c("Mean ALPS", "Left ALPS", "Right ALPS"))
)]
locus_order <- loci[order(-phenotype_count, CHR, lead_position), locus_id]
membership[, locus_label := factor(locus_label,
  levels = membership[locus_id %in% locus_order,
                      unique(locus_label), by = locus_id]$V1)]
fwrite(membership, file.path(data_dir, "Figure2C_locus_membership.csv"))

p2c <- ggplot(membership, aes(locus_label, phenotype_label, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "0" = unname(cols["verylight"]), "1" = unname(cols["green"])),
                    guide = "none") +
  labs(x = NULL, y = NULL, title = "Fourteen merged ALPS loci") +
  theme_pub(7.6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y = element_text(face = "bold")
  )

fig2 <- (p2a | p2b) / p2c +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))
save_figure(fig2, "Figure2_ALPS_genetic_architecture", 7.2, 7.0)

# ---------------------------------------------------------------------------
# Figure 3: targeted cross-trait association screen
# ---------------------------------------------------------------------------

phe <- fread(file.path(project, "04_PheWAS", "results",
                       "ALPS_targeted_PheWAS_results.tsv"), na.strings = c("", "NA"))
phe_plot <- phe[!is.na(p) & evidence_level != "Excluded"]
phe_plot[, trait_short := fifelse(trait == "Cerebral small vessel disease", "WMH",
  fifelse(trait == "Alzheimer disease", "Alzheimer\ndisease",
  fifelse(trait == "Daytime sleepiness", "Daytime\nsleepiness",
  fifelse(trait == "Sleep duration", "Sleep\nduration",
  fifelse(trait == "Type 2 diabetes", "Type 2\ndiabetes", trait)))))]
trait_order <- c("Insomnia", "Sleep\nduration", "Daytime\nsleepiness",
                 "Alzheimer\ndisease", "Parkinson disease", "Dementia",
                 "Stroke", "WMH", "Depression", "Anxiety",
                 "Body mass index", "Obesity", "Type 2\ndiabetes")
phe_plot[, trait_short := factor(trait_short, levels = trait_order)]
phe_plot[, neglog10p := -log10(p)]
phe_plot[, direction := factor(
  fifelse(beta_per_ALPS_lowering_allele >= 0, "Positive", "Negative"),
  levels = c("Positive", "Negative")
)]
phe_plot[, highlight := evidence_level %in% c("A_Bonferroni", "B_FDR")]
phe_plot[, label := fifelse(evidence_level == "A_Bonferroni", SNP, "")]
fwrite(phe_plot, file.path(data_dir, "Figure3_PheWAS_plot_data.csv"))

domain_values <- c(
  Sleep = unname(cols["blue"]), Neurodegenerative = unname(cols["purple"]),
  Cerebrovascular = unname(cols["vermilion"]),
  Psychiatric = unname(cols["orange"]), Metabolic = unname(cols["green"])
)
p3 <- ggplot(phe_plot, aes(trait_short, neglog10p)) +
  geom_hline(yintercept = -log10(0.05 / nrow(phe_plot)), linetype = 2,
             color = cols["vermilion"], linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 3,
             color = cols["grey"], linewidth = 0.45) +
  geom_jitter(
    aes(color = domain, shape = direction, alpha = highlight),
    width = 0.22, height = 0, size = 2.1, stroke = 0.35
  ) +
  ggrepel::geom_text_repel(
    data = phe_plot[label != ""],
    aes(label = label, color = domain),
    size = 2.4, max.overlaps = Inf, min.segment.length = 0,
    box.padding = 0.25, point.padding = 0.2, show.legend = FALSE
  ) +
  scale_color_manual(values = domain_values) +
  scale_shape_manual(values = c(24, 25)) +
  scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 0.95), guide = "none") +
  labs(
    x = NULL, y = expression(-log[10](italic(P))),
    color = "Trait domain", shape = "Effect of ALPS-lowering allele",
    title = "Targeted cross-trait screen of 14 merged ALPS loci"
  ) +
  theme_pub(8.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom", legend.box = "vertical"
  )
save_figure(p3, "Figure3_targeted_cross_trait_screen", 7.2, 5.6)

# ---------------------------------------------------------------------------
# Figure 4: LDSC genetic correlations
# ---------------------------------------------------------------------------

rg <- fread(file.path(project, "Gate5C", "04_LDSC_results",
                      "ALPS_LDSC_rg_results.csv"))
rg[, phenotype_label := factor(ALPS_phenotype,
  levels = c("right ALPS", "left ALPS", "mALPS"),
  labels = c("Right ALPS", "Left ALPS", "Mean ALPS")
)]
rg[, trait_label := factor(Trait,
  levels = rev(c("WMH", "Stroke", "Alzheimer disease", "Dementia", "Insomnia", "BMI"))
)]
rg[, ci_low := rg - 1.96 * SE]
rg[, ci_high := rg + 1.96 * SE]
rg[, sig_label := fifelse(Bonferroni_significant, "Bonferroni + FDR", "Not corrected-significant")]
fwrite(rg, file.path(data_dir, "Figure4_LDSC_rg_plot_data.csv"))

p4 <- ggplot(rg, aes(rg, trait_label, color = phenotype_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = cols["grey"], linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    position = position_dodge(width = 0.55), width = 0.10,
    orientation = "y", linewidth = 0.55
  ) +
  geom_point(
    aes(shape = sig_label),
    position = position_dodge(width = 0.55), size = 2.5, stroke = 0.7
  ) +
  scale_color_manual(values = c(
    "Mean ALPS" = unname(cols["blue"]), "Left ALPS" = unname(cols["orange"]),
    "Right ALPS" = unname(cols["green"])
  )) +
  scale_shape_manual(values = c("Bonferroni + FDR" = 16, "Not corrected-significant" = 1)) +
  scale_x_continuous(
    limits = c(-0.30, 0.30), breaks = seq(-0.3, 0.3, 0.1),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(
    x = expression("Genetic correlation (" * r[g] * ")"), y = NULL,
    color = NULL, shape = NULL,
    title = "Genome-wide genetic correlations with ALPS phenotypes"
  ) +
  theme_pub(9) +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "bold"))
save_figure(p4, "Figure4_LDSC_genetic_correlations", 7.2, 4.8)

# ---------------------------------------------------------------------------
# Figure 5: chr16 colocalization, SuSiE, and independent evidence
# ---------------------------------------------------------------------------

coloc <- fread(file.path(project, "Gate5B", "03_coloc_results",
                         "Gate5B_coloc_summary.csv"))
chr16 <- coloc[region_id == "ALPS_L16_011" & external_trait == "WMH"]
chr16[, phenotype_label := factor(
  fifelse(ALPS_phenotype == "mALPS", "Mean ALPS",
          fifelse(ALPS_phenotype == "left_ALPS", "Left ALPS", "Right ALPS")),
  levels = c("Mean ALPS", "Left ALPS", "Right ALPS")
)]
classic_long <- melt(
  chr16[, .(phenotype_label, PP4, PP4_strict_p12_1e_6)],
  id.vars = "phenotype_label", variable.name = "prior", value.name = "PP4_value"
)
classic_long[, prior := factor(prior,
  levels = c("PP4", "PP4_strict_p12_1e_6"),
  labels = c("p₁₂ = 10⁻⁵", "p₁₂ = 10⁻⁶")
)]
fwrite(classic_long, file.path(data_dir, "Figure5A_classic_coloc.csv"))

p5a <- ggplot(classic_long, aes(phenotype_label, PP4_value, fill = prior)) +
  geom_hline(yintercept = 0.8, linetype = 2, color = cols["vermilion"], linewidth = 0.5) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = sprintf("%.3f", PP4_value)),
            position = position_dodge(width = 0.72), vjust = -0.35, size = 2.5) +
  scale_fill_manual(values = unname(c(cols["blue"], cols["sky"]))) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "PP4", fill = "Shared prior",
       title = "Classic coloc in overlapping primary GWAS") +
  theme_pub(8) +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold"))

susie <- fread(file.path(project, "Gate5C", "05_finemapping",
                         "Gate5C_SuSiE_coloc_signal_pairs.csv"))
susie[, phenotype_label := factor(
  fifelse(ALPS_phenotype == "mALPS", "Mean ALPS",
          fifelse(ALPS_phenotype == "left_ALPS", "Left ALPS", "Right ALPS")),
  levels = c("Mean ALPS", "Left ALPS", "Right ALPS")
)]
susie_long <- melt(
  susie[, .(phenotype_label, PP.H3.abf, PP.H4.abf)],
  id.vars = "phenotype_label", variable.name = "hypothesis", value.name = "posterior"
)
susie_long[, hypothesis := factor(hypothesis,
  levels = c("PP.H3.abf", "PP.H4.abf"), labels = c("PP3", "PP4"))]
fwrite(susie_long, file.path(data_dir, "Figure5B_SuSiE_coloc.csv"))

p5b <- ggplot(susie_long, aes(phenotype_label, posterior, fill = hypothesis)) +
  geom_col(width = 0.68) +
  geom_text(
    data = susie_long[hypothesis == "PP4"],
    aes(label = sprintf("%.3f", posterior)), vjust = 1.4,
    color = "white", fontface = "bold", size = 2.6
  ) +
  scale_fill_manual(values = c(
    "PP3" = unname(cols["lightgrey"]), "PP4" = unname(cols["purple"]))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Posterior probability", fill = NULL,
       title = "Multi-signal coloc with SuSiE") +
  theme_pub(8) +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold"))

replication <- fread(file.path(project, "Gate5C", "01_WMH_independent",
                               "Gate5C_independent_WMH_coloc_summary.csv"))
replication[, phenotype_label := factor(
  fifelse(ALPS_phenotype == "mALPS", "Mean ALPS",
          fifelse(ALPS_phenotype == "left_ALPS", "Left ALPS", "Right ALPS")),
  levels = c("Mean ALPS", "Left ALPS", "Right ALPS")
)]
primary_rep <- rbind(
  chr16[, .(phenotype_label, dataset = "Primary WMH (overlapping)\nN=48,454", PP4_value = PP4)],
  replication[, .(phenotype_label, dataset = "Traylor non-UKB\nN=3,670", PP4_value = PP4)]
)
primary_rep[, dataset := factor(dataset,
  levels = c("Primary WMH (overlapping)\nN=48,454", "Traylor non-UKB\nN=3,670"))]
fwrite(primary_rep, file.path(data_dir, "Figure5C_replication_coloc.csv"))

p5c <- ggplot(primary_rep, aes(phenotype_label, PP4_value, fill = dataset)) +
  geom_hline(yintercept = 0.8, linetype = 2, color = cols["vermilion"], linewidth = 0.5) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_fill_manual(values = unname(c(cols["blue"], cols["orange"]))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(
    x = NULL, y = "PP4", fill = NULL,
    title = "Formal colocalization evidence by dataset",
    subtitle = "CHARGE regional screen: minimum P=1.22×10⁻⁶; PP4 not estimable"
  ) +
  theme_pub(8) +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold"))

cs <- fread(file.path(project, "Gate5C", "05_finemapping",
                      "Gate5C_SuSiE_credible_sets.csv"))
cs_summary <- cs[, .(
  trait, lead_SNP, lead_PIP, credible_set_size,
  credible_set_coverage = coverage
)]
cs_summary[, trait_label := factor(
  fifelse(trait == "WMH", "WMH",
          fifelse(trait == "mALPS", "Mean ALPS",
          fifelse(trait == "left_ALPS", "Left ALPS", "Right ALPS"))),
  levels = rev(c("WMH", "Mean ALPS", "Left ALPS", "Right ALPS"))
)]
fwrite(cs_summary, file.path(data_dir, "Figure5D_finemapping_summary.csv"))

p5d <- ggplot(cs_summary, aes(lead_PIP, trait_label, size = credible_set_size,
                              color = trait_label)) +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = paste0(lead_SNP, "\nCS=", credible_set_size)),
            nudge_x = 0.025, hjust = 0, size = 2.5, color = cols["black"]) +
  scale_color_manual(values = c(
    "WMH" = unname(cols["vermilion"]), "Mean ALPS" = unname(cols["blue"]),
    "Left ALPS" = unname(cols["orange"]), "Right ALPS" = unname(cols["green"])
  ), guide = "none") +
  scale_size_continuous(range = c(3, 7), guide = "none") +
  scale_x_continuous(limits = c(0, 0.36), breaks = seq(0, 0.3, 0.1)) +
  labs(x = "Lead-variant posterior inclusion probability", y = NULL,
       title = "SuSiE credible signals") +
  theme_pub(8)

fig5 <- (p5a | p5b) / (p5c | p5d) +
  plot_annotation(
    tag_levels = "A",
    caption = paste0(
      "Primary PP4: highly overlapping GWAS; internal, overlap-sensitive regional evidence.\n",
      "SuSiE models multiple signals but not sample-overlap covariance."
    ),
    theme = theme(
      plot.tag = element_text(face = "bold"),
      plot.caption = element_text(size = 7.5, color = cols["black"], hjust = 0)
    )
  )
save_figure(fig5, "Figure5_chr16_colocalization_finemapping_replication", 7.2, 7.0)

# The Q2 redesign replaces the compact posterior-only panel above with stacked
# regional association tracks, SuSiE credible-set positions, a complete p12
# sensitivity grid, and a transparent external-evaluation panel.
source(file.path(
  project, "09_manuscript", "scripts", "build_figure5_q2.R"
))

# ---------------------------------------------------------------------------
# Figure 6: chr16 functional annotation
# ---------------------------------------------------------------------------

credible <- fread(file.path(project, "Gate5D", "01_credible_variants",
                            "chr16_SuSiE_credible_variant_union.csv"))
shared_count <- credible[n_trait_credible_sets == 4, .N]

functional_counts <- data.table(
  metric = factor(
    c("Credible-set\nvariant union", "Shared by WMH +\n3 ALPS traits",
      "Mapped to\nGTEx v8", "Significant all-tissue\neQTL records",
      "Significant brain\neQTL records"),
    levels = rev(c("Credible-set\nvariant union", "Shared by WMH +\n3 ALPS traits",
                   "Mapped to\nGTEx v8", "Significant all-tissue\neQTL records",
                   "Significant brain\neQTL records"))
  ),
  count = c(nrow(credible), shared_count, 37, 59, 0),
  group = c("Variant", "Variant", "Regulatory", "Regulatory", "Regulatory")
)
fwrite(functional_counts, file.path(data_dir, "Figure6A_annotation_counts.csv"))

p6a <- ggplot(functional_counts, aes(count, metric, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = count), hjust = -0.25, size = 2.7) +
  scale_fill_manual(values = c(
    "Variant" = unname(cols["blue"]), "Regulatory" = unname(cols["orange"]))) +
  scale_x_continuous(limits = c(0, 66), expand = expansion(mult = c(0, 0))) +
  labs(x = "Count", y = NULL, fill = NULL,
       title = "Credible variants and regulatory annotation") +
  theme_pub(8) +
  theme(legend.position = "top")

genes <- fread(file.path(project, "Gate5D", "04_gene_prioritization",
                         "chr16_gene_prioritization.csv"))
gene_plot <- genes[!is.na(minimum_MAGMA_P) & gene %in%
                     c("C16orf95", "FBXO31", "ZCCHC14", "MAP1LC3B", "JPH3", "KLHDC4")]
gene_plot <- melt(
  gene_plot[, .(gene, mALPS_MAGMA_P, left_ALPS_MAGMA_P, right_ALPS_MAGMA_P)],
  id.vars = "gene", variable.name = "phenotype", value.name = "P"
)
gene_plot[, phenotype := factor(phenotype,
  levels = c("mALPS_MAGMA_P", "left_ALPS_MAGMA_P", "right_ALPS_MAGMA_P"),
  labels = c("Mean ALPS", "Left ALPS", "Right ALPS"))]
gene_plot[, gene := factor(gene, levels = rev(c(
  "C16orf95", "FBXO31", "ZCCHC14", "MAP1LC3B", "JPH3", "KLHDC4"
)))]
gene_plot[, neglog10p := -log10(P)]
fwrite(gene_plot, file.path(data_dir, "Figure6B_chr16_gene_results.csv"))

p6b <- ggplot(gene_plot, aes(neglog10p, gene, color = phenotype)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2,
             color = cols["grey"], linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.55), size = 2.4) +
  scale_color_manual(values = c(
    "Mean ALPS" = unname(cols["blue"]), "Left ALPS" = unname(cols["orange"]),
    "Right ALPS" = unname(cols["green"])
  )) +
  labs(x = expression(-log[10](italic(P))), y = NULL, color = NULL,
       title = "MAGMA results for chr16 candidate genes") +
  theme_pub(8) +
  theme(legend.position = "top", axis.text.y = element_text(face = "italic"))

pathways <- fread(file.path(project, "Gate5D", "05_pathways",
                            "ALPS_MAGMA_pathway_cross_trait_convergence.csv"))
path_plot <- pathways[traits_P_lt_0_005 == 3][1:6]
path_plot[, term_short := c(
  "CRMP/Sema3A signaling", "Microtubule regulation",
  "Proteasomal catabolism", "Calcium ion import",
  "Smooth-muscle differentiation", "Chromatin organization"
)]
path_plot[, term_short := factor(term_short, levels = rev(term_short))]
path_plot[, neglog10p := -log10(maximum_P)]
fwrite(path_plot, file.path(data_dir, "Figure6C_pathway_results.csv"))

p6c <- ggplot(path_plot, aes(neglog10p, term_short)) +
  geom_vline(xintercept = -log10(0.005), linetype = 2,
             color = cols["grey"], linewidth = 0.5) +
  geom_col(fill = cols["purple"], width = 0.65) +
  geom_text(aes(label = sprintf("%.1e", maximum_P)), hjust = -0.15, size = 2.5) +
  scale_x_continuous(limits = c(0, max(path_plot$neglog10p) + 0.8)) +
  labs(x = expression(-log[10]("maximum " * italic(P) * " across ALPS phenotypes")),
       y = NULL, title = "Nominal cross-phenotype pathway convergence") +
  theme_pub(8)

evidence <- data.table(
  gene = factor(c("C16orf95", "FBXO31", "MAP1LC3B", "ZCCHC14", "JPH3"),
                levels = rev(c("C16orf95", "FBXO31", "MAP1LC3B", "ZCCHC14", "JPH3"))),
  `Regional position` = c(1, 1, 1, 1, 1),
  `Credible-variant overlap` = c(1, 0, 0, 0, 0),
  `Brain eQTL` = c(0, 0, 0, 0, 0),
  `MAGMA nominal in 3 traits` = c(1, 0, 0, 0, 0)
)
evidence_long <- melt(evidence, id.vars = "gene",
                      variable.name = "evidence", value.name = "present")
fwrite(evidence_long, file.path(data_dir, "Figure6D_gene_evidence_matrix.csv"))

p6d <- ggplot(evidence_long, aes(evidence, gene, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_manual(values = c(
    "0" = unname(cols["verylight"]), "1" = unname(cols["green"])),
                    guide = "none") +
  labs(x = NULL, y = NULL, title = "Evidence used for regional gene prioritization") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(face = "italic")
  )

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))
save_figure(fig6, "Figure6_chr16_functional_annotation", 7.2, 7.3)

# ---------------------------------------------------------------------------
# Supplementary Figure S1: exploratory mALPS -> WMH MR
# ---------------------------------------------------------------------------

mr <- fread(file.path(project, "Gate5C", "10_mALPS_to_WMH_MR", "03_results",
                      "mALPS_to_Traylor2016_WMH_MR_results.csv"))
mr[, method_short := factor(
  fifelse(method == "Inverse variance weighted", "IVW",
  fifelse(method == "Robust adjusted profile score (RAPS)", "MR-RAPS", method)),
  levels = rev(c("IVW", "Weighted median", "MR Egger", "MR-RAPS"))
)]
fwrite(mr, file.path(data_dir, "SupplementaryFigureS1A_MR_methods.csv"))

p7a <- ggplot(mr, aes(b, method_short)) +
  geom_vline(xintercept = 0, linetype = 2, color = cols["grey"], linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.12,
                orientation = "y", color = cols["blue"], linewidth = 0.65) +
  geom_point(size = 2.8, color = cols["blue"]) +
  labs(x = "Effect on WMH burden per unit higher mALPS", y = NULL,
       title = "MR estimates") +
  theme_pub(8.5)

loo <- fread(file.path(project, "Gate5C", "10_mALPS_to_WMH_MR", "04_sensitivity",
                       "leave_one_out_results.csv"))
loo[, SNP := factor(SNP, levels = rev(SNP))]
fwrite(loo, file.path(data_dir, "SupplementaryFigureS1B_leave_one_out.csv"))

p7b <- ggplot(loo, aes(b, SNP)) +
  geom_vline(xintercept = 0, linetype = 2, color = cols["grey"], linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.08,
                orientation = "y", color = cols["green"], linewidth = 0.55) +
  geom_point(size = 2.0, color = cols["green"]) +
  labs(x = "IVW estimate after excluding each SNP", y = NULL,
       title = "Leave-one-out analysis") +
  theme_pub(7.8)

diag <- data.table(
  test = factor(
    c("IVW Cochran Q", "MR-Egger Cochran Q", "Egger intercept", "MR-PRESSO global"),
    levels = rev(c("IVW Cochran Q", "MR-Egger Cochran Q",
                   "Egger intercept", "MR-PRESSO global"))
  ),
  P = c(0.722, 0.648, 0.680, 0.762)
)
fwrite(diag, file.path(data_dir, "SupplementaryFigureS1C_diagnostics.csv"))

p7c <- ggplot(diag, aes(P, test)) +
  geom_vline(xintercept = 0.05, linetype = 2, color = cols["vermilion"], linewidth = 0.55) +
  geom_segment(aes(x = 0, xend = P, yend = test), color = cols["lightgrey"], linewidth = 1.0) +
  geom_point(size = 3, color = cols["orange"]) +
  geom_text(aes(label = sprintf("P=%.3f", P)), nudge_x = 0.025,
            hjust = 0, size = 2.6) +
  scale_x_continuous(limits = c(0, 0.98), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "P value", y = NULL, title = "Heterogeneity and pleiotropy diagnostics") +
  theme_pub(8)

strength <- data.table(
  metric = factor(c("Minimum F", "Median F", "Steiger-supporting SNPs"),
                  levels = rev(c("Minimum F", "Median F", "Steiger-supporting SNPs"))),
  value = c(29.84, 37.32, 9),
  display = c("29.84", "37.32", "9/10")
)
fwrite(strength, file.path(data_dir, "SupplementaryFigureS1D_strength_direction.csv"))

p7d <- ggplot(strength, aes(1, metric)) +
  geom_tile(width = 1.7, height = 0.70, fill = unname(cols["verylight"]),
            color = "white", linewidth = 0.8) +
  geom_text(aes(label = display), size = 4.0, fontface = "bold",
            color = unname(cols["purple"])) +
  scale_x_continuous(limits = c(0, 2), breaks = NULL) +
  labs(x = NULL, y = NULL, title = "Instrument strength and directionality") +
  theme_pub(8) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

fig7 <- (p7a | p7c) / (p7b | p7d) +
  plot_layout(widths = c(1.15, 1), heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))
save_figure(fig7, "Supplementary_Figure_S1_exploratory_mALPS_WMH_MR", 7.2, 7.2)

# ---------------------------------------------------------------------------
# Main-text tables
# ---------------------------------------------------------------------------

table1 <- data.table(
  Analysis_role = c(
    "Exposure/imaging phenotype", "Primary WMH comparison",
    "Non-UKB formal regional evaluation and supplementary MR outcome",
    "Non-UKB regional WMH screen"
  ),
  Dataset = c(
    "UK Biobank DTI-ALPS GWAS (Huang et al. 2025)",
    "Sargurupremraj et al. 2020 European WMH GWAS",
    "Traylor et al. 2016 stroke-cohort WMH GWAS",
    "Verhaaren et al. 2015 CHARGE European WMH GWAS"
  ),
  Phenotype = c(
    "Mean, left, and right DTI-ALPS indices",
    "Quantitative WMH burden",
    "Quantitative WMH burden in ischemic-stroke cohorts",
    "MRI WMH burden, ln(WMH + 1)"
  ),
  Sample_size = c("Up to 31,021", "48,454", "3,670", "Up to 17,936"),
  Ancestry = c("White British", "European", "European", "European"),
  UKB_overlap = c(
    "ALPS source cohort; up to 31,021 UKB imaging participants",
    "Includes 26,788 UKB imaging participants; maximum potential overlap 86.4%",
    "No UKB participants; full independence from non-UKB meta-analysis cohorts not established",
    "No UKB participants; full independence from non-UKB meta-analysis cohorts not established"
  ),
  Primary_use = c(
    "Genetic architecture and ALPS exposure data",
    "Targeted screen, LDSC, and overlap-sensitive primary chr16 colocalization",
    "Non-UKB chr16 colocalization evaluation and complementary directional MR evidence",
    "Regional association screen; formal PP4 not estimable"
  )
)

table2 <- merge(
  qc[, .(
    phenotype,
    GW_significant = gws_raw_rows,
    QC_retained = gws_biallelic_unique_candidates,
    Independent_sentinels = strict_independent_sentinels
  )],
  h2[, .(phenotype, h2 = h2_observed, h2_SE = SE, h2_P = P)],
  by = "phenotype"
)
rg_wmh <- rg[Trait == "WMH", .(
  phenotype = fifelse(ALPS_phenotype == "mALPS", "mALPS",
                      fifelse(ALPS_phenotype == "left ALPS", "left_ALPS", "right_ALPS")),
  WMH_rg = rg, WMH_rg_SE = SE, WMH_rg_P = P, WMH_rg_FDR = FDR_q
)]
classic_tab <- chr16[, .(
  phenotype = ALPS_phenotype,
  Classic_PP4 = PP4,
  Strict_PP4 = PP4_strict_p12_1e_6,
  Classic_lead = top_shared_candidate_SNP,
  H4_CS_size = H4_95pct_credible_set_SNPs
)]
susie_tab <- susie[, .(
  phenotype = ALPS_phenotype,
  SuSiE_PP4 = PP.H4.abf,
  SuSiE_ALPS_lead = hit1,
  SuSiE_WMH_lead = hit2
)]
table2 <- Reduce(function(x, y) merge(x, y, by = "phenotype", all = TRUE),
                 list(table2, rg_wmh, classic_tab, susie_tab))
table2[, Phenotype := fifelse(phenotype == "mALPS", "Mean ALPS",
                       fifelse(phenotype == "left_ALPS", "Left ALPS", "Right ALPS"))]
table2[, phenotype_order := match(
  Phenotype, c("Mean ALPS", "Left ALPS", "Right ALPS")
)]
setorder(table2, phenotype_order)
setcolorder(table2, c(
  "Phenotype", "GW_significant", "QC_retained", "Independent_sentinels",
  "h2", "h2_SE", "h2_P", "WMH_rg", "WMH_rg_SE", "WMH_rg_P", "WMH_rg_FDR",
  "Classic_PP4", "Strict_PP4", "Classic_lead", "H4_CS_size",
  "SuSiE_PP4", "SuSiE_ALPS_lead", "SuSiE_WMH_lead", "phenotype",
  "phenotype_order"
))
table2[, c("phenotype", "phenotype_order") := NULL]

table3 <- mr[, .(
  Method = fifelse(method == "Inverse variance weighted", "IVW",
            fifelse(method == "Robust adjusted profile score (RAPS)", "MR-RAPS", method)),
  SNPs = nsnp,
  Beta = b,
  SE = se,
  CI_lower = ci_lower,
  CI_upper = ci_upper,
  P = pval
)]
table3_diagnostics <- data.table(
  Diagnostic = c(
    "Minimum F statistic", "Median F statistic", "IVW Cochran Q",
    "MR-Egger Cochran Q", "MR-Egger intercept", "MR-PRESSO global test",
    "MR-PRESSO outliers", "Steiger-supporting variants"
  ),
  Estimate = c("29.84", "37.32", "6.180 (df=9)", "5.997 (df=8)",
               "0.0195 (SE=0.0454)", "RSSobs=7.377", "0", "9/10"),
  P = c(NA, NA, 0.722, 0.648, 0.680, 0.762, NA, NA)
)

fwrite(table1, file.path(table_dir, "Table1_data_sources.csv"))
fwrite(table2, file.path(table_dir, "Table2_core_genetic_results.csv"))
# MR estimates are supplementary and are assembled in Supplementary Table 12.

wb <- createWorkbook()
header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
body_style <- createStyle(valign = "top", wrapText = TRUE)
for (item in list(
  list("Table1_Data_sources", table1),
  list("Table2_Core_results", table2)
)) {
  addWorksheet(wb, item[[1]])
  writeData(wb, item[[1]], item[[2]], headerStyle = header_style)
  addStyle(wb, item[[1]], body_style, rows = 2:(nrow(item[[2]]) + 1),
           cols = 1:ncol(item[[2]]), gridExpand = TRUE)
  freezePane(wb, item[[1]], firstRow = TRUE)
  setColWidths(wb, item[[1]], cols = 1:ncol(item[[2]]), widths = "auto")
  setColWidths(wb, item[[1]], cols = 1:ncol(item[[2]]), widths = pmin(
    35, pmax(10, sapply(item[[2]], function(z) max(nchar(as.character(z)), na.rm = TRUE) + 2))
  ))
  addFilter(wb, item[[1]], rows = 1, cols = 1:ncol(item[[2]]))
}
saveWorkbook(wb, file.path(table_dir, "DTI_ALPS_WMH_main_text_tables.xlsx"),
             overwrite = TRUE)

cat("Created figures and tables in:", asset_dir, "\n")
