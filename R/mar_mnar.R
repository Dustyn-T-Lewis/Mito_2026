# =============================================================================
# mar_mnar.R  --  3-method MAR/MNAR consensus classification (reporting only)
#
# Votes MNAR via: (1) k-means on (mean intensity, % missing); (2) global logistic
# P(missing | intensity); (3) left-tail proximity (fraction of observed values
# below the global 25th percentile, weighted by % missing). >=2/3 votes = MNAR.
# DEP uses the NON-imputed matrix; this classification drives the imputation_reliable
# flag + reporting. Refs: Lazar 2016 (PMID 26906401), Karpievitch 2012.
#
# NOTE: kmeans is stochastic -- the CALLER must set.seed() for reproducibility.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

classify_mar_mnar <- function(mat, miss_unreliable = 50) {
  n_samp    <- ncol(mat)
  prot_miss <- rowSums(is.na(mat))
  pct_miss  <- prot_miss / n_samp * 100
  obs_means <- rowMeans(mat, na.rm = TRUE)
  has_na    <- which(prot_miss > 0 & prot_miss < n_samp)
  inc_mean  <- obs_means[has_na]
  inc_pct   <- pct_miss[has_na]

  # (1) k-means: the lower-intensity cluster votes MNAR
  km      <- kmeans(scale(cbind(inc_mean, inc_pct)), centers = 2, nstart = 25)
  km_mnar <- km$cluster == which.min(tapply(inc_mean, km$cluster, mean))

  # (2) global logistic P(missing | intensity)
  lr_fit  <- glm(is_miss ~ intensity, family = binomial,
                 data = data.frame(is_miss = as.integer(is.na(as.vector(mat))),
                                   intensity = rep(obs_means, n_samp)))
  lr_pred <- predict(lr_fit, newdata = data.frame(intensity = inc_mean), type = "response")
  lr_mnar <- lr_pred > median(lr_pred)

  # (3) left-tail proximity
  q25       <- quantile(mat, 0.25, na.rm = TRUE)
  tail_frac <- vapply(has_na, function(i) mean(mat[i, !is.na(mat[i, ])] < q25), numeric(1))
  lt_mnar   <- (tail_frac * inc_pct / 100) > median(tail_frac * inc_pct / 100)

  votes <- as.integer(km_mnar) + as.integer(lr_mnar) + as.integer(lt_mnar)

  out <- tibble(
    uniprot_id = rownames(mat), n_miss = prot_miss, pct_miss = pct_miss,
    mean_intensity = obs_means,
    vote_kmeans = NA_integer_, vote_logistic = NA_integer_,
    vote_lefttail = NA_integer_, n_mnar_votes = NA_integer_
  )
  out$vote_kmeans[has_na]   <- as.integer(km_mnar)
  out$vote_logistic[has_na] <- as.integer(lr_mnar)
  out$vote_lefttail[has_na] <- as.integer(lt_mnar)
  out$n_mnar_votes[has_na]  <- votes

  out |> mutate(
    classification = case_when(
      n_miss == 0        ~ "Complete",
      n_miss >= n_samp   ~ "MNAR",
      n_mnar_votes >= 2  ~ "MNAR",
      TRUE               ~ "MAR"
    ),
    imputation_reliable = classification == "Complete" | pct_miss < miss_unreliable
  )
}
