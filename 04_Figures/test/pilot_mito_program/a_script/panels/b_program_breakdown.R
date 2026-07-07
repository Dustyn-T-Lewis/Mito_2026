# Panel B: which mitochondrial programs change. Split MitoCarta by MitoPathways
# top level and separate structural machinery (OXPHOS, mtDNA/translation) from
# regulatory programs (metabolism, dynamics/mitophagy, import, transport,
# signaling). Bulk cargo lifts structure uniformly; a regulated response also
# moves the regulatory arms.

CAT_LABELS <- c(
  OXPHOS = "OXPHOS",
  MITOCHONDRIAL_CENTRAL_DOGMA = "mtDNA & translation",
  METABOLISM = "Metabolism",
  MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE = "Dynamics & mitophagy",
  PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS = "Import & proteostasis",
  SIGNALING = "Signaling",
  SMALL_MOLECULE_TRANSPORT = "Small-molecule transport"
)
STRUCTURAL_CATS <- c("OXPHOS", "MITOCHONDRIAL_CENTRAL_DOGMA")

build_program_breakdown <- function(comb, cat_map) {
  d <- comb |>
    dplyr::transmute(gene = toupper(gene), Transplant = logFC_CTLvMITO, Rescue = logFC_PHEvPHE_MITO) |>
    dplyr::inner_join(cat_map, by = "gene") |>
    tidyr::pivot_longer(c(Transplant, Rescue), names_to = "contrast", values_to = "logFC") |>
    dplyr::group_by(category, contrast) |>
    dplyr::summarise(median_logFC = median(logFC, na.rm = TRUE), n = dplyr::n_distinct(gene), .groups = "drop") |>
    dplyr::mutate(
      arm = ifelse(category %in% STRUCTURAL_CATS, "structural", "regulatory"),
      label = dplyr::coalesce(CAT_LABELS[category], category)
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(median_logFC, reorder(label, median_logFC), color = contrast)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::facet_grid(arm ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = c(Transplant = "#4393C3", Rescue = "#4DAF4A"), name = NULL) +
    ggplot2::labs(
      title = "Mito program breakdown",
      subtitle = "structural cargo vs regulatory programs (median over each program)",
      x = "median log2FC", y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(legend.position = "bottom")

  list(plot = p, table = dplyr::arrange(d, arm, category, contrast))
}
