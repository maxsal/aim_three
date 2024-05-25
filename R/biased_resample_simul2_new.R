ms::libri(
  ms, qs, data.table, tidyverse, glue, cli, MASS,
  survey, mice, patchwork, furrr, parallelly
)


if (availableCores() <= 8) {
    plan(multisession, workers = 6)
} else {
    wkrs <- min(c(availableCores() / 2, 18))
    plan(multicore, workers = wkrs)
}

source("fn/mat_and_means.R")
source("fn/regression_helpers.R")
source("fn/resample_helpers.R")

outcome      <- "glucose"
exposure     <- "bmi"
covs         <- c("age", "female", "nhw", "smoke")
sample_sizes <- c(1000, 2500, 5000, 10000)
iterations   <- 1000
f            <- paste0(outcome, " ~ ", exposure)
f_cov        <- paste0(f, " + ", paste0(covs, collapse = " + "))

formals(map_dfr)$.progress   <- formals(future_map)$.progress <- formals(future_map_dfr)$.progress <- formals(map)$.progress <- TRUE
formals(future_map)$.options <- formals(future_map_dfr)$.options <- furrr_options(seed = TRUE)
formals(mice)$printFlag      <- FALSE

for (sample_size in sample_sizes) {
  cli_progress_step(paste0("sample size: ", sample_size))
  res <- future_map_dfr(
    seq_len(iterations),
    \(i) {
      set.seed(i)
      data <- data.table::copy(generate_data(n = 100000, mu = means, mat = mat))

      # sampling
      samp_score <- data[, bmi] + data[, age] + data[, glucose]
      samp_prob <- plogis(samp_score)
      data[, ww := samp_prob]
      data[, wgt := 1 / ww]
      data$wgt <- data$wgt * (nrow(data) / sum(data$wgt))
      samp_idx <- sample(seq_len(nrow(data)), sample_size, prob = data[, ww])
      data$select <- 0
      data$select[samp_idx] <- 1

      # samples
      mcar_data <- data[select == 1, ]
      mcar_ex_ind <- rbinom(sample_size, c(0,1), 0.25)
      mcar_out_ind <- rbinom(sample_size, c(0,1), 0.25)
      mcar_data[mcar_ex_ind == 1, bmi := NA]
      mcar_data[mcar_out_ind == 1, glucose := NA]
      mcar_data$wgt <- mcar_data$wgt * (nrow(mcar_data) / sum(mcar_data$wgt))

      mar_data <- make_ex_out_mar3(
                  data       = data[select == 1, ],
                  exposure = exposure,
                  outcome = outcome,
                  covs = covs,
                  ex_int = -4.07,
                  out_int = -4.07
      )
      mar_data$wgt <- mar_data$wgt * (nrow(mar_data) / sum(mar_data$wgt))

      mnar_data <- make_ex_out_mnar3(
          data = data[select == 1, ],
          exposure = exposure,
          outcome = outcome,
          covs = covs,
          ex_int = -5.7,
          out_int = -5.7
      )
      mnar_data$wgt <- mnar_data$wgt * (nrow(mnar_data) / sum(mnar_data$wgt))

      # cli_progress_step("Imputing data...")
      # initialize predictor matrix
      these_vars <- c(outcome, exposure, covs)
      not_these_vars <- names(data)[!names(data) %in% these_vars]
      pred_mat <- suppressWarnings(mice(mcar_data, m = 1, printFlag = FALSE)$predictorMatrix)
      pred_mat[not_these_vars, ] <- 0
      pred_mat[, not_these_vars] <- 0

      mcar_imp <- mice(data = mcar_data, predictorMatrix = pred_mat)
      mar_imp <- mice(data = mar_data, predictorMatrix = pred_mat)
      mnar_imp <- mice(data = mnar_data, predictorMatrix = pred_mat)

      mcar_dsn <- svydesign(ids = ~1, data = mcar_data, weights = ~wgt)
      mar_dsn <- svydesign(ids = ~1, data = mar_data, weights = ~wgt)
      mnar_dsn <- svydesign(ids = ~1, data = mnar_data, weights = ~wgt)

      # cli_progress_step("Imputing data with PRS...")
      these_vars_prs <- c(these_vars, "bmi_prs", "glucose_prs")
      not_these_vars_prs <- names(data)[!names(data) %in% these_vars_prs]
      pred_mat_prs <- suppressWarnings(mice(mcar_data, m = 1, printFlag = FALSE)$predictorMatrix)
      pred_mat_prs[not_these_vars_prs, ] <- 0
      pred_mat_prs[, not_these_vars_prs] <- 0

      mcar_imp_prs <- mice(data = mcar_data, predictorMatrix = pred_mat_prs)
      mar_imp_prs <- mice(data = mar_data, predictorMatrix = pred_mat_prs)
      mnar_imp_prs <- mice(data = mnar_data, predictorMatrix = pred_mat_prs)

      rbindlist(list(
        # truth
        .glm(as.formula(f), data = data, var = exposure)[, `:=`(method = "Full", adj = "Unadjusted", miss_data = NA_character_, weight = NA_character_)],
        .glm(as.formula(f_cov), data = data, var = exposure)[, `:=`(method = "Full", adj = "Adjusted", miss_data = NA_character_, weight = NA_character_)],

        # mcar
        ## unadjusted
        .glm(as.formula(f), data = mcar_data, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MCAR", weight = "Unweighted")],
        .glm_pool(as.formula(f), mcar_imp, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MCAR", weight = "Unweighted")],
        .glm_pool(as.formula(f), mcar_imp_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MCAR", weight = "Unweighted")],
        .svyglm(as.formula(f), design = mcar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MCAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f), imp_object = mcar_imp, weights = mcar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MCAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f), imp_object = mcar_imp_prs, weights = mcar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MCAR", weight = "Weighted")],
        ## adjusted
        .glm(as.formula(f_cov), data = mcar_data, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MCAR", weight = "Unweighted")],
        .glm_pool(as.formula(f_cov), mcar_imp, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MCAR", weight = "Unweighted")],
        .glm_pool(as.formula(f_cov), mcar_imp_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MCAR", weight = "Unweighted")],
        .svyglm(as.formula(f_cov), design = mcar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MCAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f_cov), imp_object = mcar_imp, weights = mcar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MCAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f_cov), imp_object = mcar_imp_prs, weights = mcar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MCAR", weight = "Weighted")],

        # mar
        ## unadjusted
        .glm(as.formula(f), data = mar_data, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MAR", weight = "Unweighted")],
        .glm_pool(as.formula(f), mar_imp, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MAR", weight = "Unweighted")],
        .glm_pool(as.formula(f), mar_imp_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MAR", weight = "Unweighted")],
        .svyglm(as.formula(f), design = mar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f), imp_object = mar_imp, weights = mar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f), imp_object = mar_imp_prs, weights = mar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MAR", weight = "Weighted")],
        ## adjusted
        .glm(as.formula(f_cov), data = mar_data, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MAR", weight = "Unweighted")],
        .glm_pool(as.formula(f_cov), mar_imp, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MAR", weight = "Unweighted")],
        .glm_pool(as.formula(f_cov), mar_imp_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MAR", weight = "Unweighted")],
        .svyglm(as.formula(f_cov), design = mar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f_cov), imp_object = mar_imp, weights = mar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f_cov), imp_object = mar_imp_prs, weights = mar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MAR", weight = "Weighted")],

        # mnar
        ## unadjusted
        .glm(as.formula(f), data = mnar_data, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MNAR", weight = "Unweighted")],
        .glm_pool(as.formula(f), mnar_imp, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MNAR", weight = "Unweighted")],
        .glm_pool(as.formula(f), mnar_imp_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MNAR", weight = "Unweighted")],
        .svyglm(as.formula(f), design = mnar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MNAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f), imp_object = mnar_imp, weights = mnar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MNAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f), imp_object = mnar_imp_prs, weights = mnar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MNAR", weight = "Weighted")],
        ## adjusted
        .glm(as.formula(f_cov), data = mnar_data, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MNAR", weight = "Unweighted")],
        .glm_pool(as.formula(f_cov), mnar_imp, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MNAR", weight = "Unweighted")],
        .glm_pool(as.formula(f_cov), mnar_imp_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MNAR", weight = "Unweighted")],
        .svyglm(as.formula(f_cov), design = mnar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MNAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f_cov), imp_object = mnar_imp, weights = mnar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MNAR", weight = "Weighted")],
        .svyglm_pool(as.formula(f_cov), imp_object = mnar_imp_prs, weights = mnar_data[["wgt"]], var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MNAR", weight = "Weighted")]
      ), use.names = TRUE, fill = TRUE)[, `:=`(seed = i, size = sample_size)][]
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  # res <- rbindlist(res, use.names = TRUE, fill = TRUE)
  fwrite(
    x    = res,
    file = paste0("data/public/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_biased_simul_res.csv")
  )

  # cli_progress_step("calculating diagnostics...")
  truth_una <- res[method == "Full" & adj == "Unadjusted", ]
  truth_adj <- res[method == "Full" & adj == "Adjusted", ]

  mcar_cc_unw_una <- res[method == "Complete case" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mcar_imp_unw_una <- res[method == "Imputed" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mcar_impp_unw_una <- res[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mcar_cc_w_una <- res[method == "Complete case" & miss_data == "MCAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mcar_imp_w_una <- res[method == "Imputed" & miss_data == "MCAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mcar_impp_w_una <- res[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Weighted" & adj == "Unadjusted", ]

  mcar_cc_unw_adj <- res[method == "Complete case" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mcar_imp_unw_adj <- res[method == "Imputed" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mcar_impp_unw_adj <- res[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mcar_cc_w_adj <- res[method == "Complete case" & miss_data == "MCAR" & weight == "Weighted" & adj == "Adjusted", ]
  mcar_imp_w_adj <- res[method == "Imputed" & miss_data == "MCAR" & weight == "Weighted" & adj == "Adjusted", ]
  mcar_impp_w_adj <- res[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Weighted" & adj == "Adjusted", ]

  mar_cc_unw_una <- res[method == "Complete case" & miss_data == "MAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mar_imp_unw_una <- res[method == "Imputed" & miss_data == "MAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mar_impp_unw_una <- res[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mar_cc_w_una <- res[method == "Complete case" & miss_data == "MAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mar_imp_w_una <- res[method == "Imputed" & miss_data == "MAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mar_impp_w_una <- res[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Weighted" & adj == "Unadjusted", ]

  mar_cc_unw_adj <- res[method == "Complete case" & miss_data == "MAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mar_imp_unw_adj <- res[method == "Imputed" & miss_data == "MAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mar_impp_unw_adj <- res[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mar_cc_w_adj <- res[method == "Complete case" & miss_data == "MAR" & weight == "Weighted" & adj == "Adjusted", ]
  mar_imp_w_adj <- res[method == "Imputed" & miss_data == "MAR" & weight == "Weighted" & adj == "Adjusted", ]
  mar_impp_w_adj <- res[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Weighted" & adj == "Adjusted", ]

  mnar_cc_unw_una <- res[method == "Complete case" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mnar_imp_unw_una <- res[method == "Imputed" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mnar_impp_unw_una <- res[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mnar_cc_w_una <- res[method == "Complete case" & miss_data == "MNAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mnar_imp_w_una <- res[method == "Imputed" & miss_data == "MNAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mnar_impp_w_una <- res[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Weighted" & adj == "Unadjusted", ]

  mnar_cc_unw_adj <- res[method == "Complete case" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mnar_imp_unw_adj <- res[method == "Imputed" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mnar_impp_unw_adj <- res[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mnar_cc_w_adj <- res[method == "Complete case" & miss_data == "MNAR" & weight == "Weighted" & adj == "Adjusted", ]
  mnar_imp_w_adj <- res[method == "Imputed" & miss_data == "MNAR" & weight == "Weighted" & adj == "Adjusted", ]
  mnar_impp_w_adj <- res[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Weighted" & adj == "Adjusted", ]

  diag <- rbindlist(list(
    # mcar
    ## unadjusted
    simul_diag(cleaned = mcar_cc_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mcar_imp_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mcar_impp_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mcar_cc_w_una, truth = truth_una$est),
    simul_diag(cleaned = mcar_imp_w_una, truth = truth_una$est),
    simul_diag(cleaned = mcar_impp_w_una, truth = truth_una$est),
    ## adjusted
    simul_diag(cleaned = mcar_cc_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mcar_imp_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mcar_impp_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mcar_cc_w_adj, truth = truth_adj$est),
    simul_diag(cleaned = mcar_imp_w_adj, truth = truth_adj$est),
    simul_diag(cleaned = mcar_impp_w_adj, truth = truth_adj$est),

    # mar
    ## unadjusted
    simul_diag(cleaned = mar_cc_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mar_imp_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mar_impp_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mar_cc_w_una, truth = truth_una$est),
    simul_diag(cleaned = mar_imp_w_una, truth = truth_una$est),
    simul_diag(cleaned = mar_impp_w_una, truth = truth_una$est),
    ## adjusted
    simul_diag(cleaned = mar_cc_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mar_imp_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mar_impp_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mar_cc_w_adj, truth = truth_adj$est),
    simul_diag(cleaned = mar_imp_w_adj, truth = truth_adj$est),
    simul_diag(cleaned = mar_impp_w_adj, truth = truth_adj$est),

    # mnar
    ## unadjusted
    simul_diag(cleaned = mnar_cc_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mnar_imp_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mnar_impp_unw_una, truth = truth_una$est),
    simul_diag(cleaned = mnar_cc_w_una, truth = truth_una$est),
    simul_diag(cleaned = mnar_imp_w_una, truth = truth_una$est),
    simul_diag(cleaned = mnar_impp_w_una, truth = truth_una$est),
    ## adjusted
    simul_diag(cleaned = mnar_cc_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mnar_imp_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mnar_impp_unw_adj, truth = truth_adj$est),
    simul_diag(cleaned = mnar_cc_w_adj, truth = truth_adj$est),
    simul_diag(cleaned = mnar_imp_w_adj, truth = truth_adj$est),
    simul_diag(cleaned = mnar_impp_w_adj, truth = truth_adj$est)
  ))

  fwrite(
    x = diag,
    file = paste0("data/public/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_biased_simul_diag.csv")
  )
}
