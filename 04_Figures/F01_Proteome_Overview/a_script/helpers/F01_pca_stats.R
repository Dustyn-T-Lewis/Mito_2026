# Group-structure stats for the F01 PCA: overall PERMANOVA, the pre-specified
# pairwise contrasts (BH-adjusted), and a betadisper dispersion check. All on the
# scaled Euclidean distance and seeded so the permutation p-values reproduce.
pca_group_stats <- function(mat, groups, pairs, seed = 42, permutations = 9999) {
  dist_mat <- dist(scale(t(mat)))
  set.seed(seed)
  overall <- vegan::adonis2(dist_mat ~ Group,
    data = data.frame(Group = groups),
    permutations = permutations, by = "terms"
  )

  pair_res <- dplyr::bind_rows(lapply(names(pairs), function(nm) {
    keep <- groups %in% pairs[[nm]]
    smeta <- data.frame(Group = droplevels(groups[keep]))
    set.seed(seed)
    r <- vegan::adonis2(dist(scale(t(mat[, keep]))) ~ Group,
      data = smeta, permutations = permutations, by = "terms"
    )
    tibble::tibble(role = nm, R2 = r$R2[1], p = r$`Pr(>F)`[1])
  }))
  pair_res$q <- p.adjust(pair_res$p, "BH")

  set.seed(seed)
  bd_p <- vegan::permutest(vegan::betadisper(dist_mat, groups),
    permutations = 999
  )$tab$`Pr(>F)`[1]
  if (!is.na(bd_p) && bd_p < 0.05) {
    warning("Heterogeneous group dispersions (betadisper p < 0.05)")
  }

  list(
    overall_R2 = overall["Group", "R2"],
    overall_p = overall["Group", "Pr(>F)"],
    pairwise = pair_res,
    betadisper_p = bd_p,
    table = dplyr::bind_rows(
      tibble::tibble(role = "Group (overall)", R2 = overall["Group", "R2"], p = overall["Group", "Pr(>F)"]),
      pair_res,
      tibble::tibble(role = "dispersion (betadisper)", R2 = NA_real_, p = bd_p)
    )
  )
}
