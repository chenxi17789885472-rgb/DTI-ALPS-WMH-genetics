#!/usr/bin/env Rscript

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5D"))) return(normalizePath(candidate, mustWork = TRUE))
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

root <- file.path(resolve_project_root(), "Gate5D")
region <- read.csv(file.path(root, "04_gene_prioritization", "chr16_region_MAGMA_gene_results.csv"),
                   check.names = FALSE)
conv <- read.csv(file.path(root, "05_pathways", "ALPS_MAGMA_pathway_cross_trait_convergence.csv"),
                 check.names = FALSE)

gene_order <- c("C16orf95", "FBXO31", "ZCCHC14", "MAP1LC3B", "JPH3", "KLHDC4")
trait_order <- c("mALPS", "left_ALPS", "right_ALPS")
gene_mat <- sapply(trait_order, function(tt) {
  x <- region[region$TRAIT == tt, ]
  -log10(x$P[match(gene_order, x$SYMBOL)])
})
rownames(gene_mat) <- gene_order

top <- conv[conv$traits_P_lt_0_005 == 3, ]
top <- top[order(top$maximum_P), ][1:min(8, nrow(top)), ]
short_term <- gsub("Regulation Of ", "Reg. ", top$term)
short_term <- gsub("Proteasomal Protein Catabolic Process", "Proteasomal protein catabolism", short_term)
short_term <- gsub("Smooth Muscle Cell Differentiation", "Smooth muscle cell differentiation", short_term)
short_term <- gsub("Chromatin Organization", "Chromatin organization", short_term)
short_term <- gsub("Calcium Ion Import", "Calcium ion import", short_term)

draw_plot <- function() {
  par(mfrow = c(1, 2), mar = c(8, 4.5, 4, 1), oma = c(1, 1, 3, 1),
      family = "Helvetica")
  bp <- barplot(t(gene_mat), beside = TRUE,
                col = c("#1F4E79", "#4F81BD", "#9DC3E6"),
                names.arg = gene_order, las = 2,
                ylab = expression(-log[10]("MAGMA gene P")),
                main = "chr16 regional gene evidence",
                border = NA)
  abline(h = -log10(0.05), col = "#A61C00", lty = 2, lwd = 1.2)
  legend("topright", legend = c("mALPS", "left ALPS", "right ALPS", "nominal P=0.05"),
         fill = c("#1F4E79", "#4F81BD", "#9DC3E6", NA),
         border = c(NA, NA, NA, NA), lty = c(NA, NA, NA, 2),
         col = c(NA, NA, NA, "#A61C00"), bty = "n", cex = 0.8)

  par(mar = c(5, 12, 4, 1))
  vals <- -log10(rev(top$maximum_P))
  labs <- rev(short_term)
  barplot(vals, horiz = TRUE, names.arg = labs, las = 1,
          col = "#70AD47", border = NA,
          xlab = expression(-log[10]("maximum P across three traits")),
          main = "Cross-trait nominal pathway convergence",
          cex.names = 0.72)
  abline(v = -log10(0.005), col = "#A61C00", lty = 2, lwd = 1.2)
  mtext("Exploratory: no pathway survived BH-FDR or Bonferroni correction.",
        side = 1, line = 3.6, cex = 0.72, col = "#7F6000")
  mtext("Gate 5D chr16 functional annotation", outer = TRUE, side = 3,
        line = 1, cex = 1.25, font = 2)
}

png(file.path(root, "06_figures", "Gate5D_chr16_functional_annotation.png"),
    width = 3000, height = 1450, res = 220)
draw_plot()
dev.off()

pdf(file.path(root, "06_figures", "Gate5D_chr16_functional_annotation.pdf"),
    width = 13.6, height = 6.6)
draw_plot()
dev.off()
