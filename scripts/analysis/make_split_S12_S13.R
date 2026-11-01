suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  library(grid)
})

source_path <- Sys.getenv("MR_LOO_SOURCE")
out_dir <- Sys.getenv("SUPP_FIG_OUT_DIR")
package_figure_dir <- Sys.getenv("SUPP_FIG_PACKAGE_DIR")
stopifnot(
  nzchar(source_path),
  nzchar(out_dir),
  nzchar(package_figure_dir)
)

loo <- read.csv(
  source_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM",
  colClasses = "character"
)
loo$b <- as.numeric(loo$b)
loo$se <- as.numeric(loo$se)
loo$p <- as.numeric(loo$p)
loo$or <- exp(loo$b)
loo$or_low <- exp(loo$b - 1.96 * loo$se)
loo$or_high <- exp(loo$b + 1.96 * loo$se)
loo$nominal <- loo$p < 0.05
loo$is_all <- loo$SNP == "All"

fmt_p <- function(p) {
  ifelse(
    p < 0.001,
    formatC(p, format = "e", digits = 1),
    sprintf("%.3f", p)
  )
}

make_panel <- function(label, panel, display) {
  dat <- loo[loo$outcome_label == label, , drop = FALSE]
  all_row <- dat[dat$SNP == "All", , drop = FALSE]
  snp_rows <- dat[dat$SNP != "All", , drop = FALSE]
  snp_rows <- snp_rows[order(snp_rows$b), , drop = FALSE]
  dat <- rbind(all_row, snp_rows)
  dat$row_label <- ifelse(
    dat$SNP == "All",
    "All SNPs (IVW)",
    paste0("Omit ", dat$SNP)
  )
  dat$row_label <- factor(
    dat$row_label,
    levels = rev(c("All SNPs (IVW)", paste0("Omit ", snp_rows$SNP)))
  )
  dat$value_label <- sprintf(
    "%.2f (%.2f-%.2f)\nP=%s",
    dat$or,
    dat$or_low,
    dat$or_high,
    fmt_p(dat$p)
  )

  ggplot(dat, aes(y = row_label, x = or)) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "#777777",
      linewidth = 0.4
    ) +
    geom_errorbar(
      aes(xmin = or_low, xmax = or_high, color = nominal),
      orientation = "y",
      width = 0.18,
      linewidth = 0.55
    ) +
    geom_point(aes(color = nominal, shape = is_all), size = 1.8) +
    geom_text(
      aes(x = 1.72, label = value_label),
      hjust = 0,
      size = 2.45,
      lineheight = 0.88,
      color = "#333333"
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 18),
      guide = "none"
    ) +
    scale_color_manual(
      values = c(`FALSE` = "#235789", `TRUE` = "#C84C09"),
      guide = "none"
    ) +
    scale_x_log10(
      limits = c(0.5, 6),
      breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2),
      labels = function(x) formatC(x, format = "fg", digits = 3)
    ) +
    labs(
      title = paste0(panel, " | ", display),
      x = "OR (log scale)",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 8.8, base_family = "sans") +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 9.4,
        color = "#111111"
      ),
      axis.title.x = element_text(size = 8.2, margin = margin(t = 4)),
      axis.text.x = element_text(size = 7.8),
      axis.text.y = element_text(size = 7.8),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = unit(3.0, "pt"),
      panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#F0F0F0", linewidth = 0.22),
      plot.margin = margin(8, 12, 8, 10)
    )
}

save_group <- function(
  panels,
  stem,
  heading,
  show_method_icons = FALSE,
  show_color_legend = FALSE
) {
  blank <- nullGrob()
  if (length(panels) == 3) {
    top_row <- arrangeGrob(
      panels[[1]],
      panels[[2]],
      ncol = 2,
      widths = c(1, 1)
    )
    bottom_row <- arrangeGrob(
      blank,
      panels[[3]],
      blank,
      ncol = 3,
      widths = c(1, 2, 1)
    )
    panel_grid <- arrangeGrob(
      top_row,
      bottom_row,
      ncol = 1,
      heights = c(1, 1)
    )
  } else {
    panel_grid <- arrangeGrob(
      grobs = panels,
      ncol = 2,
      nrow = 2
    )
  }
  heading_grob <- textGrob(
    heading,
    gp = gpar(
      fontsize = 12,
      fontface = "bold",
      col = "#111111"
    )
  )
  if (show_method_icons) {
    key_grob <- function(graphic, label) {
      arrangeGrob(
        graphic,
        textGrob(
          label,
          gp = gpar(fontsize = 7.5, col = "#333333")
        ),
        ncol = 2,
        widths = c(0.28, 1)
      )
    }
    method_grob <- arrangeGrob(
      key_grob(
        pointsGrob(
          0.5, 0.5, pch = 18,
          size = unit(2.4, "mm"),
          gp = gpar(col = "#235789")
        ),
        "Complete IVW fit"
      ),
      key_grob(
        pointsGrob(
          0.5, 0.5, pch = 16,
          size = unit(2.1, "mm"),
          gp = gpar(col = "#235789")
        ),
        "After SNP omission"
      ),
      key_grob(
        segmentsGrob(
          x0 = c(0.20, 0.20, 0.80),
          x1 = c(0.80, 0.20, 0.80),
          y0 = c(0.50, 0.30, 0.30),
          y1 = c(0.50, 0.70, 0.70),
          gp = gpar(col = "#235789", lwd = 1.4)
        ),
        "95% CI"
      ),
      key_grob(
        segmentsGrob(
          0.10, 0.50, 0.90, 0.50,
          gp = gpar(col = "#777777", lwd = 1.2, lty = 2)
        ),
        "OR = 1 reference"
      ),
      nrow = 1,
      widths = c(1.45, 1.55, 0.95, 1.35)
    )
    note_grob <- method_grob
  } else {
    note_parts <- c(
      "Diamond: complete IVW fit",
      "Circle: IVW after SNP omission",
      "Error bars: 95% CI",
      "Dashed line: OR = 1"
    )
    if (!show_color_legend) {
      note_parts <- c(note_parts, "Orange = nominal P < 0.05")
    }
    note_grob <- textGrob(
      paste(note_parts, collapse = "   |   "),
      gp = gpar(fontsize = 7.5, col = "#333333")
    )
  }

  if (show_color_legend) {
    legend_source <- ggplot(
      data.frame(
        status = factor(
          c("P < 0.05", "P >= 0.05"),
          levels = c("P < 0.05", "P >= 0.05")
        ),
        x = 1:2,
        y = 1
      ),
      aes(x, y, color = status)
    ) +
      geom_point(size = 2.5) +
      scale_color_manual(
        values = c("P < 0.05" = "#C84C09", "P >= 0.05" = "#235789"),
        breaks = c("P < 0.05", "P >= 0.05"),
        labels = c("Nominal P < 0.05", "Nominal P >= 0.05"),
        name = NULL
      ) +
      theme_void() +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 7.5, color = "#333333"),
        legend.key.size = unit(0.35, "cm")
      )
    legend_grob <- cowplot::get_legend(legend_source)
    if (show_method_icons) {
      header <- arrangeGrob(
        heading_grob,
        note_grob,
        legend_grob,
        ncol = 1,
        heights = c(1.0, 0.9, 0.45)
      )
    } else {
      header <- arrangeGrob(
        heading_grob,
        note_grob,
        legend_grob,
        ncol = 1,
        heights = c(1.0, 0.8, 0.45)
      )
    }
  } else {
    if (show_method_icons) {
      header <- arrangeGrob(
        heading_grob,
        note_grob,
        ncol = 1,
        heights = c(1.0, 0.9)
      )
    } else {
      header <- arrangeGrob(
        heading_grob,
        note_grob,
        ncol = 1,
        heights = c(1.1, 0.8)
      )
    }
  }
  figure_height <- if (show_color_legend) {
    9.7
  } else if (show_method_icons) {
    9.55
  } else {
    9.3
  }
  figure <- arrangeGrob(
    header,
    panel_grid,
    ncol = 1,
    heights = c(
      if (show_color_legend) 1.15 else if (show_method_icons) 1.0 else 0.9,
      9.2
    )
  )
  tmp_dir <- tempfile(paste0(stem, "_"))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(tmp_dir, paste0(stem, ".pdf"))
  png_path <- file.path(tmp_dir, paste0(stem, ".png"))
  svg_path <- file.path(tmp_dir, paste0(stem, ".svg"))
  cairo_pdf(pdf_path, width = 10, height = figure_height, onefile = TRUE)
  grid.newpage()
  grid.draw(figure)
  dev.off()
  png(
    png_path,
    width = 10 * 600,
    height = figure_height * 600,
    res = 600,
    bg = "white",
    type = "cairo"
  )
  grid.newpage()
  grid.draw(figure)
  dev.off()
  svg(svg_path, width = 10, height = figure_height, onefile = TRUE)
  grid.newpage()
  grid.draw(figure)
  dev.off()
  stopifnot(
    file.copy(pdf_path, file.path(out_dir, paste0(stem, ".pdf")), overwrite = TRUE),
    file.copy(png_path, file.path(out_dir, paste0(stem, ".png")), overwrite = TRUE),
    file.copy(svg_path, file.path(out_dir, paste0(stem, ".svg")), overwrite = TRUE)
  )
}

figure_s12 <- list(
  make_panel("FinnGen_NAFLD", "A", "FinnGen NAFLD"),
  make_panel("FinnGen_CIRRHOSIS_BROAD", "B", "FinnGen cirrhosis"),
  make_panel(
    "FinnGen_K11_FIBROCHIRLIV",
    "C",
    "FinnGen fibrotic liver disease"
  ),
  make_panel(
    "GCST90054782_NAFLD_4761cases",
    "D",
    "NAFLD (4,761 cases)"
  )
)
figure_s13 <- list(
  make_panel(
    "GCST90091033_NAFLD_8434cases",
    "A",
    "NAFLD (8,434 cases)"
  ),
  make_panel(
    "GCST90319877_CIRRHOSIS_15225cases",
    "B",
    "Cirrhosis (15,225 cases)"
  ),
  make_panel(
    "GCST90728570_MASLD_MTAG_N122644",
    "C",
    "MASLD MTAG (N = 122,644)"
  )
)

save_group(
  figure_s12,
  "FigureS12",
  "Supplementary Figure S12 | Leave-one-out sensitivity (part 1)",
  show_method_icons = TRUE,
  show_color_legend = TRUE
)
save_group(
  figure_s13,
  "FigureS13",
  "Supplementary Figure S13 | Leave-one-out sensitivity (part 2)",
  show_method_icons = TRUE,
  show_color_legend = FALSE
)

file.copy(
  file.path(out_dir, "FigureS12.png"),
  file.path(package_figure_dir, "FigureS12_leave_one_out_part1.png"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, "FigureS12.pdf"),
  file.path(package_figure_dir, "FigureS12_leave_one_out_part1.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, "FigureS12.svg"),
  file.path(package_figure_dir, "FigureS12_leave_one_out_part1.svg"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, "FigureS13.png"),
  file.path(package_figure_dir, "FigureS13_leave_one_out_part2.png"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, "FigureS13.pdf"),
  file.path(package_figure_dir, "FigureS13_leave_one_out_part2.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, "FigureS13.svg"),
  file.path(package_figure_dir, "FigureS13_leave_one_out_part2.svg"),
  overwrite = TRUE
)

cat("Generated split Supplementary Figures S12 and S13\n")
