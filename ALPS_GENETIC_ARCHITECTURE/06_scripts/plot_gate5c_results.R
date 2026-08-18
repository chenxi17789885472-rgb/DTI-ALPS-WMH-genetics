suppressPackageStartupMessages(library(data.table))

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5C"))) return(normalizePath(candidate, mustWork = TRUE))
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

root <- resolve_project_root()
fig_dir <- file.path(root, "Gate5C/06_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rg <- fread(file.path(root, "Gate5C/04_LDSC_results/ALPS_LDSC_rg_results.csv"))
trait_order <- c("WMH", "Stroke", "Alzheimer disease", "Dementia", "Insomnia", "BMI")
pheno_order <- c("mALPS", "left ALPS", "right ALPS")
rg[, y := match(Trait, rev(trait_order))]
offsets <- c("mALPS" = 0.22, "left ALPS" = 0, "right ALPS" = -0.22)
cols <- c("mALPS" = "#2457A7", "left ALPS" = "#D17A00", "right ALPS" = "#2E8B57")

draw_forest <- function() {
  par(mar = c(4.5, 8.5, 2.3, 1.5), las = 1)
  plot(
    NA, xlim = c(-0.34, 0.24), ylim = c(0.5, 6.5),
    xlab = "Genetic correlation (rg)", ylab = "", yaxt = "n",
    main = "Genome-wide genetic correlations with DTI-ALPS"
  )
  abline(v = 0, lty = 2, col = "grey55")
  axis(2, at = 1:6, labels = rev(trait_order), tick = FALSE)
  for (p in pheno_order) {
    d <- rg[ALPS_phenotype == p]
    yy <- d$y + offsets[[p]]
    segments(d$rg - 1.96 * d$SE, yy, d$rg + 1.96 * d$SE, yy, col = cols[[p]], lwd = 2)
    points(d$rg, yy, pch = ifelse(d$Bonferroni_significant, 17, 16),
           col = cols[[p]], cex = 1.05)
  }
  legend("bottomright", legend = pheno_order, col = cols[pheno_order],
         pch = 19, bty = "n", cex = 0.85)
}

png(file.path(fig_dir, "Gate5C_LDSC_forest.png"), width = 1800, height = 1150, res = 180)
draw_forest()
dev.off()
pdf(file.path(fig_dir, "Gate5C_LDSC_forest.pdf"), width = 9, height = 5.8)
draw_forest()
dev.off()

primary <- fread(file.path(root, "Gate5B/03_coloc_results/Gate5B_coloc_summary.csv"))
primary <- primary[external_trait == "WMH", .(ALPS_phenotype, PP4_primary_ABF = PP4)]
ind <- fread(file.path(root, "Gate5C/01_WMH_independent/Gate5C_independent_WMH_coloc_summary.csv"))
ind <- ind[, .(ALPS_phenotype, PP4_independent_ABF = PP4)]
susie <- fread(file.path(root, "Gate5C/05_finemapping/Gate5C_SuSiE_coloc_signal_pairs.csv"))
susie <- susie[, .(ALPS_phenotype, PP4_SuSiE = PP.H4.abf)]
rob <- Reduce(function(x, y) merge(x, y, by = "ALPS_phenotype"), list(primary, ind, susie))
rob <- rob[match(c("mALPS", "left_ALPS", "right_ALPS"), ALPS_phenotype)]
mat <- t(as.matrix(rob[, .(PP4_primary_ABF, PP4_independent_ABF, PP4_SuSiE)]))
colnames(mat) <- c("mALPS", "left ALPS", "right ALPS")
rownames(mat) <- c("Combined WMH: ABF", "Non-UKB WMH: ABF", "Combined WMH: SuSiE")

draw_coloc <- function() {
  par(mar = c(4.5, 5, 2.3, 1.2))
  barplot(
    mat, beside = TRUE, ylim = c(0, 1.08),
    col = c("#2457A7", "#B8C2CC", "#2E8B57"),
    ylab = "Posterior probability of shared signal (PP4)",
    main = "chr16 ALPS–WMH colocalization robustness",
    legend.text = rownames(mat),
    args.legend = list(x = "topright", bty = "n", cex = 0.78)
  )
  abline(h = 0.8, lty = 2, col = "#A33A2B", lwd = 1.5)
}
png(file.path(fig_dir, "Gate5C_coloc_robustness.png"), width = 1700, height = 1100, res = 180)
draw_coloc()
dev.off()
pdf(file.path(fig_dir, "Gate5C_coloc_robustness.pdf"), width = 8.5, height = 5.5)
draw_coloc()
dev.off()

pip <- fread(cmd = paste("gzip -dc", shQuote(file.path(root, "Gate5C/05_finemapping/Gate5C_SuSiE_variant_PIP.csv.gz"))))
draw_pip <- function() {
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.4, 0.8))
  for (p in c("mALPS", "left_ALPS", "right_ALPS")) {
    d <- pip[ALPS_phenotype == p]
    plot(
      d$PIP_ALPS, d$PIP_WMH, pch = 16, cex = 0.45,
      col = adjustcolor("#2457A7", alpha.f = 0.38),
      xlab = paste0(p, " PIP"), ylab = "WMH PIP",
      main = p, xlim = c(0, max(d$PIP_ALPS) * 1.08),
      ylim = c(0, max(d$PIP_WMH) * 1.08)
    )
    top <- d[order(-(PIP_ALPS * PIP_WMH))][1:3]
    points(top$PIP_ALPS, top$PIP_WMH, pch = 19, col = "#A33A2B")
    text(top$PIP_ALPS, top$PIP_WMH, labels = top$SNP, pos = 4, cex = 0.65)
  }
}
png(file.path(fig_dir, "Gate5C_SuSiE_PIP.png"), width = 2200, height = 800, res = 180)
draw_pip()
dev.off()
pdf(file.path(fig_dir, "Gate5C_SuSiE_PIP.pdf"), width = 11, height = 4)
draw_pip()
dev.off()
