ms::libri(
  ms, qs, data.table, tidyverse, glue, cli, MASS,
  survey, mice, patchwork, furrr, parallelly
)

means <- c(
  glucose = -1.96825433778349e-16, bmi = -7.91723213769746e-17,
  age = 1.20111356353818e-16, female = 1.42619696156583e-16, nhw = -5.03312408648097e-18,
  smoke = -3.2111110600884e-17, bmi_prs = -1.15519592499013e-14,
  glucose_prs = -7.07320752607373e-15
)
# correlation matrix
# mat <- structure(
#     c(
#         1, -0.138472771069173, 0.262539473913268, 0.0468405451550121,
#         0.167566551860689, 0.166092379443715, -0.0556810150535702, -0.00389196515634526,
#         -0.138472771069173, 1, -0.179053565661372, 0.0235312040482812,
#         -0.115388276481859, -0.0448208957108339, 0.0237946498434344,
#         -0.00768675768547915, 0.262539473913268, -0.179053565661372,
#         1, 0.203021425445865, 0.0915425306467291, 0.0580461239663158,
#         0.0826686083769788, 0.0905772810900326, 0.0468405451550121, 0.0235312040482812,
#         0.203021425445865, 1, 0.0204487766742965, -0.00921594984143435,
#         0.303774286296125, -0.00405822627168262, 0.167566551860689, -0.115388276481859,
#         0.0915425306467291, 0.0204487766742965, 1, 0.0801551548165317,
#         0.0743051790132458, 0.00156663252186455, 0.166092379443715, -0.0448208957108339,
#         0.0580461239663158, -0.00921594984143435, 0.0801551548165317,
#         1, -0.0365680107481083, 0.0227178773800331, -0.0556810150535702,
#         0.0237946498434344, 0.0826686083769788, 0.303774286296125, 0.0743051790132458,
#         -0.0365680107481083, 1, 0.0636833004109739, -0.00389196515634526,
#         -0.00768675768547915, 0.0905772810900326, -0.00405822627168262,
#         0.00156663252186455, 0.0227178773800331, 0.0636833004109739,
#         1
#     ),
#     dim = c(8L, 8L),
#     dimnames = list(
#         c(
#             "age", "female", "glucose", "bmi", "smoke", "nhw",
#             "bmi_prs", "glucose_prs"
#         ),
#         c(
#             "age", "female", "glucose", "bmi", "smoke", "nhw",
#             "bmi_prs", "glucose_prs"
#         )
#     )
# ) |> as.matrix

# covariance matrix
mat <- structure(c(
  0.997523893105382, 0.208644976915785, 0.253708389862574,
  -0.177297309127456, 0.0596898614979812, 0.0917005397091668, 0.0838985048101617,
  0.0922531746254993, 0.208644976915785, 0.996951208183969, 0.0506946067208541,
  0.0203494958260491, -0.0055661377527309, 0.0216563774116052,
  0.300569363860713, -0.00500820212770039, 0.253708389862574, 0.0506946067208541,
  0.923159608174341, -0.126399831360388, 0.160140678440441, 0.165362429492189,
  -0.0519840639263659, -0.00358944704781315, -0.177297309127456,
  0.0203494958260491, -0.126399831360388, 1.00138532594553, -0.0381787297773077,
  -0.114761906189592, 0.0220751526556959, -0.00805369111191262,
  0.0596898614979812, -0.0055661377527309, 0.160140678440441, -0.0381787297773077,
  0.90677551228785, 0.0775607612957751, -0.0328217604862329, 0.0226931865599238,
  0.0917005397091668, 0.0216563774116052, 0.165362429492189, -0.114761906189592,
  0.0775607612957751, 1.00513020712039, 0.075917841346984, 0.00115950734428159,
  0.0838985048101617, 0.300569363860713, -0.0519840639263659, 0.0220751526556959,
  -0.0328217604862329, 0.075917841346984, 0.990270282102947, 0.0612015767047263,
  0.0922531746254993, -0.00500820212770039, -0.00358944704781315,
  -0.00805369111191262, 0.0226931865599238, 0.00115950734428159,
  0.0612015767047263, 1.00474034266984
), dim = c(8L, 8L), dimnames = list(
  c(
    "glucose", "bmi", "age", "female", "nhw", "smoke", "bmi_prs",
    "glucose_prs"
  ), c(
    "glucose", "bmi", "age", "female", "nhw",
    "smoke", "bmi_prs", "glucose_prs"
  )
)) |> as.matrix()

select_threshold <- function(data, stable, variable, step = 0.05) {
  tmp_seq <- seq(min(as.numeric(data[[variable]])) + step, max(as.numeric(data[[variable]])) - step, step)
  cors <- map_dbl(
    tmp_seq,
    \(i) {
      cor(data[[stable]], as.numeric(data[[variable]] >= i))
    }
  )
  cor_diff <- cors - cor(data[[stable]], data[[variable]])
  tmp_seq[which.min(abs(cor_diff))]
}

generate_data <- function(
  n        = 100000,
  mu       = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
  bin_vars = c("female", "smoke", "nhw"),
  mat      = NULL
) {
  data <- MASS::mvrnorm(n = n,
        mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
        Sigma = mat) |> data.table::as.data.table()
  data[, id := 1:.N]
  if (!is.null(bin_vars)) {
    for (bin_var in bin_vars) {
      data[[bin_var]] <- as.numeric(data[[bin_var]] >= select_threshold(data, "age", bin_var))
    }
  }
  data
}

map_generate_data <- function(
    iterations = 1000,
    size = 100000,
    mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
    mat = NULL) {
  furrr::future_map(
    seq_len(iterations),
    \(x) {
      set.seed(x)
      generate_data(n = size, mu = mu, mat = mat)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
}

  # simul_diag <- function(cleaned, truth) {
  #   raw_bias <- mean(cleaned$est - truth)
  #   per_bias <- mean(100 * abs(raw_bias / truth))
  #   cov_rate <- mean(cleaned$lo < truth & cleaned$hi > truth)
  #   avg_wide <- mean(cleaned$hi - cleaned$lo)
  #   rmse <- sqrt(mean((cleaned$est - truth)^2))
  #   data.table(
  #     adj       = cleaned[, unique(adj)],
  #     miss_data = cleaned[, unique(miss_data)],
  #     weight    = cleaned[, unique(weight)],
  #     size      = cleaned[, unique(size)],
  #     method    = cleaned[, unique(method)],
  #     raw_bias  = raw_bias,
  #     per_bias  = per_bias,
  #     cov_rate  = cov_rate,
  #     avg_wide  = avg_wide,
  #     rmse      = rmse
  #   )
  # }

if (availableCores() <= 8) {
  plan(multisession, workers = 6)
} else {
  wkrs <- min(c(availableCores() / 2, 24))
  plan(multicore, workers = wkrs)
}


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
      mcar_ex_ind <- rbinom(nrow(mcar_data), 1, 0.25)
      mcar_out_ind <- rbinom(nrow(mcar_data), 1, 0.25)
      mcar_data[mcar_ex_ind == 1, exposure] <- NA
      mcar_data[mcar_out_ind == 1, outcome] <- NA
      mcar_data$wgt <- mcar_data$wgt * (nrow(mcar_data) / sum(mcar_data$wgt))

      mar_data <- make_ex_out_mar5(
        data = data[select == 1, ],
        exposure = exposure,
        outcome = outcome,
        covs = covs,
        ex_int = -3.73,
        out_int = -3.13
      )
      mar_data$wgt <- mar_data$wgt * (nrow(mar_data) / sum(mar_data$wgt))

      mnar_data <- make_ex_out_mnar5(
          data = data[select == 1, ],
          exposure = exposure,
          outcome = outcome,
          covs = covs,
          ex_int = -4.52,
          out_int = -3.2
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
    file = paste0("data/public/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_exoutmiss_biased_simul_res.csv")
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
    file = paste0("data/public/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_exoutmiss_biased_simul_diag.csv")
  )
}
