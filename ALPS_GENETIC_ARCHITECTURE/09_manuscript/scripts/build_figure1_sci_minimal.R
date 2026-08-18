#!/usr/bin/env Rscript

# Figure contract
# Core conclusion: DTI-ALPS and WMH are MRI-derived phenotypes evaluated through
# a summary-statistics imaging-genetics workflow.
# Evidence chain: A, DTI-ALPS measurement concept; B, real FLAIR WMH example;
# C, analysis sequence from GWAS to shared architecture and triangulation.
# Archetype: image plate + schematic-led composite.
# Backend: R only.
# Export: 183 x 125 mm; editable SVG/PDF; 600 dpi TIFF; 300 dpi PNG preview.

# Rscript does not reliably expose the script path through commandArgs()[1].
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
script_path <- normalizePath(file_arg, mustWork = TRUE)
manuscript_dir <- dirname(dirname(script_path))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
  library(svglite)
  library(ragg)
})

source_png <- file.path(
  manuscript_dir,
  "imaging_figure_sources",
  "Cai_2022_Figure1_FLAIR_PWMH_DWMH_CC_BY_crop.png"
)
figure_dir <- file.path(manuscript_dir, "submission_assets", "figures")
archive_dir <- file.path(figure_dir, "archive")
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

base_name <- file.path(figure_dir, "Figure1_study_design")
old_files <- paste0(base_name, c(".svg", ".pdf", ".png", ".tiff"))

# Preserve the immediately preceding version once, without changing its pixels.
for (ext in c("svg", "pdf", "png", "tiff")) {
  old_file <- paste0(base_name, ".", ext)
  archive_file <- file.path(archive_dir, paste0("Figure1_ornate_v3.", ext))
  if (file.exists(old_file) && !file.exists(archive_file)) {
    file.copy(old_file, archive_file, overwrite = FALSE)
  }
}

font_family <- "Helvetica"
ink <- "#222222"
mid_grey <- "#666666"
light_grey <- "#D9D9D9"
very_light <- "#F3F3F3"
accent <- "#2F6F9F"
accent_light <- "#EAF1F6"

theme_figure <- theme_void(base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

# Panel A: original, restrained DTI-ALPS schematic.
theta <- seq(0, 2 * pi, length.out = 300)
brain <- data.frame(
  x = 0.35 + 0.25 * cos(theta),
  y = 0.51 + 0.34 * sin(theta)
)

# Build ventricular silhouettes as simple polygons.
vent_l <- data.frame(
  x = c(0.305, 0.282, 0.278, 0.292, 0.318, 0.332, 0.326),
  y = c(0.62, 0.57, 0.49, 0.43, 0.45, 0.53, 0.60)
)
vent_r <- transform(vent_l, x = 0.70 - x)

p_a <- ggplot() +
  geom_path(data = brain, aes(x, y), linewidth = 0.45, colour = mid_grey) +
  geom_polygon(data = vent_l, aes(x, y), fill = "#4D4D4D", colour = NA) +
  geom_polygon(data = vent_r, aes(x, y), fill = "#4D4D4D", colour = NA) +
  # Projection fibers: solid vertical lines
  annotate(
    "segment",
    x = c(0.22, 0.25, 0.45, 0.48),
    xend = c(0.22, 0.25, 0.45, 0.48),
    y = c(0.30, 0.27, 0.27, 0.30),
    yend = c(0.72, 0.75, 0.75, 0.72),
    linewidth = 0.75,
    colour = "#777777",
    lineend = "round"
  ) +
  # Association fibers: dashed arcs approximated by polylines
  annotate(
    "path",
    x = c(0.18, 0.14, 0.13, 0.15, 0.19),
    y = c(0.69, 0.60, 0.49, 0.39, 0.31),
    linewidth = 0.55,
    linetype = "22",
    colour = mid_grey
  ) +
  annotate(
    "path",
    x = c(0.52, 0.56, 0.57, 0.55, 0.51),
    y = c(0.69, 0.60, 0.49, 0.39, 0.31),
    linewidth = 0.55,
    linetype = "22",
    colour = mid_grey
  ) +
  # Perivascular direction and four ROIs: two projection and two association ROIs.
  annotate(
    "segment",
    x = c(0.175, 0.39),
    xend = c(0.285, 0.50),
    y = c(0.51, 0.51),
    yend = c(0.51, 0.51),
    linewidth = 0.75,
    colour = accent
  ) +
  annotate(
    "point",
    x = c(0.23, 0.47),
    y = c(0.51, 0.51),
    shape = 21,
    size = 3.8,
    stroke = 0.8,
    fill = "white",
    colour = accent
  ) +
  annotate(
    "point",
    x = c(0.15, 0.55),
    y = c(0.51, 0.51),
    shape = 22,
    size = 3.6,
    stroke = 0.8,
    fill = "white",
    colour = mid_grey
  ) +
  # Direct labels
  annotate(
    "segment",
    x = 0.66,
    xend = 0.72,
    y = 0.69,
    yend = 0.69,
    colour = "#777777",
    linewidth = 0.75
  ) +
  annotate(
    "text",
    x = 0.74,
    y = 0.69,
    label = "Projection fibers (z)",
    hjust = 0,
    size = 6.5 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "point",
    x = 0.69,
    y = 0.38,
    shape = 21,
    size = 2.7,
    stroke = 0.65,
    fill = "white",
    colour = accent
  ) +
  annotate(
    "text",
    x = 0.74,
    y = 0.38,
    label = "Projection ROIs (n = 2)",
    hjust = 0,
    size = 6.5 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "point",
    x = 0.69,
    y = 0.29,
    shape = 22,
    size = 2.6,
    stroke = 0.65,
    fill = "white",
    colour = mid_grey
  ) +
  annotate(
    "text",
    x = 0.74,
    y = 0.29,
    label = "Association ROIs (n = 2)",
    hjust = 0,
    size = 6.5 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "segment",
    x = 0.66,
    xend = 0.72,
    y = 0.59,
    yend = 0.59,
    colour = mid_grey,
    linetype = "22",
    linewidth = 0.55
  ) +
  annotate(
    "text",
    x = 0.74,
    y = 0.59,
    label = "Association fibers (y)",
    hjust = 0,
    size = 6.5 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "segment",
    x = 0.66,
    xend = 0.72,
    y = 0.49,
    yend = 0.49,
    colour = accent,
    linewidth = 0.75
  ) +
  annotate(
    "text",
    x = 0.74,
    y = 0.49,
    label = "Perivascular direction (x)",
    hjust = 0,
    size = 6.5 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 0.92,
    label = "A",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 8 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.09,
    y = 0.92,
    label = "DTI-ALPS measurement",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 7 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.50,
    y = -0.025,
    label = "ALPS index = mean(Dxproj, Dxassoc) / mean(Dyproj, Dzassoc)",
    hjust = 0.5,
    size = 6.1 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.35,
    y = 0.085,
    label = "Four spherical ROIs (radius = 4 mm) at the lateral-ventricular level",
    hjust = 0.5,
    size = 5.8 / .pt,
    family = font_family,
    colour = mid_grey
  ) +
  coord_fixed(xlim = c(0, 1.15), ylim = c(-0.08, 0.95), clip = "off") +
  theme_figure

# Panel B: real FLAIR/PWMH/DWMH image, with no decorative container.
mri <- png::readPNG(source_png)
mri_grob <- rasterGrob(mri, interpolate = TRUE)

p_b <- ggplot() +
  annotation_custom(mri_grob, xmin = 0.03, xmax = 0.97, ymin = 0.12, ymax = 0.78) +
  annotate(
    "text",
    x = 0.02,
    y = 0.94,
    label = "B",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 8 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.10,
    y = 0.94,
    label = "WMH on axial T2-FLAIR MRI",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 7 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.50,
    y = 0.04,
    label = "Illustrative external image; analyses used total WMH burden",
    hjust = 0.5,
    size = 5.8 / .pt,
    family = font_family,
    colour = mid_grey
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_figure

# Panel C: condensed linear workflow, with a single highlighted core-analysis step.
workflow <- data.frame(
  x = 1:6,
  label = c(
    "ALPS GWAS",
    "Loci + SNP h²",
    "Targeted cross-trait\nscreen",
    "LDSC genetic\ncorrelation",
    "chr16 coloc\n+ SuSiE",
    "Non-UKB evaluation\n+ annotation"
  )
)

connections <- data.frame(x = 1:5, xend = 2:6)

p_c <- ggplot() +
  geom_segment(
    data = connections,
    aes(x = x + 0.35, xend = xend - 0.35, y = 0.48, yend = 0.48),
    linewidth = 0.45,
    colour = mid_grey,
    arrow = arrow(length = unit(1.7, "mm"), type = "closed")
  ) +
  geom_rect(
    data = workflow,
    aes(
      xmin = x - 0.34,
      xmax = x + 0.34,
      ymin = 0.28,
      ymax = 0.68,
      fill = x == 5,
      colour = x == 5
    ),
    linewidth = 0.45
  ) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = accent_light), guide = "none") +
  scale_colour_manual(values = c(`FALSE` = "#8A8A8A", `TRUE` = accent), guide = "none") +
  geom_text(
    data = workflow,
    aes(x = x, y = 0.48, label = label),
    size = 6.2 / .pt,
    lineheight = 0.95,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.56,
    y = 0.98,
    label = "C",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 8 / .pt,
    family = font_family,
    colour = ink
  ) +
  annotate(
    "text",
    x = 0.82,
    y = 0.98,
    label = "Summary-statistics imaging-genetics workflow",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 7 / .pt,
    family = font_family,
    colour = ink
  ) +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0.1, 1), expand = FALSE, clip = "off") +
  theme_figure +
  theme(plot.margin = margin(1, 3, 1, 3, unit = "mm"))

figure <- (p_a | p_b) / p_c +
  plot_layout(widths = c(1.08, 0.92), heights = c(2.0, 0.78)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

width_mm <- 183
height_mm <- 125
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(
  paste0(base_name, ".svg"),
  width = width_in,
  height = height_in,
  system_fonts = list(sans = font_family)
)
print(figure)
dev.off()

grDevices::cairo_pdf(
  paste0(base_name, ".pdf"),
  width = width_in,
  height = height_in,
  family = font_family,
  bg = "white"
)
print(figure)
dev.off()

ragg::agg_tiff(
  paste0(base_name, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(figure)
dev.off()

ragg::agg_png(
  paste0(base_name, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(figure)
dev.off()

# R-only grayscale QA preview.
preview <- png::readPNG(paste0(base_name, ".png"))
rgb <- preview[, , 1:3, drop = FALSE]
grey <- 0.2126 * rgb[, , 1] + 0.7152 * rgb[, , 2] + 0.0722 * rgb[, , 3]
grey_rgb <- array(rep(grey, 3), dim = c(dim(grey), 3))
png::writePNG(grey_rgb, paste0(base_name, "_grayscale_QA.png"))

cat("Saved:", paste0(base_name, c(".svg", ".pdf", ".tiff", ".png")), sep = "\n")
cat("\n")
