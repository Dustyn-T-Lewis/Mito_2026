# Supplementary S1: network construction diagnostics. The scale-free fit justifies the
# soft power (lowest power clearing R^2 >= 0.85), shown beside mean connectivity; the gene
# dendrogram with the merged module colours shows the clustering the modules come from.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

panel_scale_free <- function(sft_df, chosen_power, rsq_cut = 0.85) {
  d <- sft_df |>
    transmute(power = Power, rsq = SFT.R.sq, mean_k = mean.k.)
  fit <- ggplot(d, aes(power, rsq)) +
    geom_hline(yintercept = rsq_cut, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = chosen_power, color = "#B2182B", linewidth = 0.3) +
    geom_line(color = "grey55", linewidth = 0.3) +
    geom_text(aes(label = power), size = 1.8, fontface = "bold", color = "grey20") +
    labs(title = "Scale-free fit", x = "soft power", y = expression("scale-free R"^2)) +
    FIG_THEME
  conn <- ggplot(d, aes(power, mean_k)) +
    geom_vline(xintercept = chosen_power, color = "#B2182B", linewidth = 0.3) +
    geom_line(color = "grey55", linewidth = 0.3) +
    geom_text(aes(label = power), size = 1.8, fontface = "bold", color = "grey20") +
    labs(title = "Mean connectivity", x = "soft power", y = "mean k") +
    FIG_THEME
  fit | conn
}

panel_dendro <- function(net) {
  block <- net$blockGenes[[1]]
  ggplotify::as.ggplot(function() {
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]], net$colors[block],
      groupLabels = "module", dendroLabels = FALSE, addGuide = TRUE,
      hang = 0.03, guideHang = 0.05, main = "Protein dendrogram and modules",
      cex.colorLabels = 0.6, cex.main = 0.8, cex.axis = 0.6, cex.lab = 0.6, marAll = c(1, 4, 2, 0)
    )
  }) +
    theme(plot.margin = margin(t = 6, r = 2, b = 2, l = 2, unit = "mm"))
}
