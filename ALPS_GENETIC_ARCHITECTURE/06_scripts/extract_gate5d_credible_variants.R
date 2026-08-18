suppressPackageStartupMessages(library(data.table))

resolve_project_root <- function() {
  env_root <- Sys.getenv("ALPS_PROJECT_ROOT", unset = "")
  candidates <- c(env_root, file.path(getwd(), "ALPS_GENETIC_ARCHITECTURE"), getwd())
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "Gate5C")) && dir.exists(file.path(candidate, "Gate5D"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate ALPS_GENETIC_ARCHITECTURE. Run from the workspace/project root or set ALPS_PROJECT_ROOT.")
}

root <- resolve_project_root()
fm <- file.path(root, "Gate5C/05_finemapping")
out_dir <- file.path(root, "Gate5D/01_credible_variants")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("WMH", "mALPS", "left_ALPS", "right_ALPS")
files <- c(
  WMH = "susie_WMH.rds",
  mALPS = "susie_mALPS.rds",
  left_ALPS = "susie_left_ALPS.rds",
  right_ALPS = "susie_right_ALPS.rds"
)

rows <- list()
for (trait in traits) {
  fit <- readRDS(file.path(fm, files[[trait]]))
  cs <- fit$sets$cs
  if (length(cs) == 0) next
  for (i in seq_along(cs)) {
    idx <- cs[[i]]
    rows[[length(rows) + 1]] <- data.table(
      trait = trait,
      signal = names(cs)[i],
      SNP = colnames(fit$alpha)[idx],
      PIP = fit$pip[idx],
      credible_set_coverage = fit$sets$coverage[i]
    )
  }
}

long <- rbindlist(rows)
wide <- dcast(long, SNP ~ trait, value.var = "PIP", fill = 0)
wide[, max_PIP := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = traits]
wide[, n_trait_credible_sets := rowSums(.SD > 0), .SDcols = traits]
setorder(wide, -n_trait_credible_sets, -max_PIP)

fwrite(long, file.path(out_dir, "chr16_SuSiE_credible_set_membership_long.csv"))
fwrite(wide, file.path(out_dir, "chr16_SuSiE_credible_variant_union.csv"))
fwrite(wide[, .(SNP)], file.path(out_dir, "chr16_SuSiE_credible_variant_rsids.txt"), col.names = FALSE)

cat("Union credible variants:", nrow(wide), "\n")
cat("Shared by all four traits:", sum(wide$n_trait_credible_sets == 4), "\n")
print(head(wide, 20))
