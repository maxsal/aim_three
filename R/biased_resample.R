ms::libri(
  ms, qs, data.table, tidyverse, glue, cli, MASS,
  survey, mice, patchwork, furrr, parallelly
)

cor_structure <- structure(
  c(1, -0.138472771069173, 0.262539473913268, 0.0468405451550121, 
    0.167566551860689, 0.166092379443715, -0.0556810150535702, -0.00389196515634526, 
    -0.138472771069173, 1, -0.179053565661372, 0.0235312040482812, 
    -0.115388276481859, -0.0448208957108339, 0.0237946498434344, 
    -0.00768675768547915, 0.262539473913268, -0.179053565661372, 
    1, 0.203021425445865, 0.0915425306467291, 0.0580461239663158, 
    0.0826686083769788, 0.0905772810900326, 0.0468405451550121, 0.0235312040482812, 
    0.203021425445865, 1, 0.0204487766742965, -0.00921594984143435, 
    0.303774286296125, -0.00405822627168262, 0.167566551860689, -0.115388276481859, 
    0.0915425306467291, 0.0204487766742965, 1, 0.0801551548165317, 
    0.0743051790132458, 0.00156663252186455, 0.166092379443715, -0.0448208957108339, 
    0.0580461239663158, -0.00921594984143435, 0.0801551548165317, 
    1, -0.0365680107481083, 0.0227178773800331, -0.0556810150535702, 
    0.0237946498434344, 0.0826686083769788, 0.303774286296125, 0.0743051790132458, 
    -0.0365680107481083, 1, 0.0636833004109739, -0.00389196515634526, 
    -0.00768675768547915, 0.0905772810900326, -0.00405822627168262, 
    0.00156663252186455, 0.0227178773800331, 0.0636833004109739, 
    1),
    dim = c(8L, 8L),
    dimnames = list(c("age", "female", "glucose", "bmi", "smoke", "nhw",
                      "bmi_prs", "glucose_prs"),
                    c("age", "female", "glucose", "bmi", "smoke", "nhw",
                      "bmi_prs", "glucose_prs")
                    )
)

mat <- cor_structure |> as.matrix()

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

if (availableCores() <= 8) {
  plan(multisession, workers = 6)
} else if (availableCores() >= 64) {
  plan(multicore, workers = 18)
}

source("fn/regression_helpers.R")
source("fn/resample_helpers.R")

outcome      <- "glucose"
exposure     <- "bmi"
covs         <- c("age", "female", "nhw", "smoke")
sample_sizes <- c(1000, 2500, 5000, 10000)
iterations   <- 1000

f     <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(f, " + ", paste0(covs, collapse = " + "))

for (sample_size in sample_sizes) {
  cli_progress_step(glue("sample size: {sample_size}"))
  tmp <- future_map_dfr(
    1:iterations, \(i) {
      set.seed(i) 

      # generate data
      data <- generate_data(n = 100000, mat = mat)

      # selection bias and missingness
      ## selection
      score <- data$bmi + data$age + data$glucose
      prob  <- plogis(score)

      data$ww  <- prob
      data$wgt <- 1 / prob
      data$wgt <- data$wgt * (nrow(data) / sum(data$wgt))

      samp_idx              <- sample(seq_len(nrow(data)), sample_size, prob = data$ww)
      data$select           <- 0
      data$select[samp_idx] <- 1

      ## missingness
      mcar_data <- make_missing(
        data       = data,
        target_var = exposure,
        mech       = "MCAR"
      )
      mar_data <- make_missing(
        data       = data,
        target_var = exposure,
        vars       = c("age", outcome),
        var_wgts   = c(2, 2),
        mech       = "MAR",
        scale_cont = FALSE
      )
      mnar_data <- make_missing(
        data       = data,
        target_var = exposure,
        vars       = c("age", outcome),
        var_wgts   = c(2, 2),
        mech       = "MNAR",
        scale_cont = FALSE
      )

      # select sample
      mcar_data_sub <- mcar_data[samp_idx, ]
      mar_data_sub  <- mar_data[samp_idx, ]
      mnar_data_sub <- mnar_data[samp_idx, ]

      # scale weights
      mcar_data_sub$wgt <- mcar_data_sub$wgt * (nrow(mcar_data_sub) / sum(mcar_data_sub$wgt))
      mcar_data_sub_wgt <- mcar_data_sub$wgt

      mar_data_sub$wgt <- mar_data_sub$wgt * (nrow(mar_data_sub) / sum(mar_data_sub$wgt))
      mar_data_sub_wgt <- mar_data_sub$wgt

      mnar_data_sub$wgt <- mnar_data_sub$wgt * (nrow(mnar_data_sub) / sum(mnar_data_sub$wgt))
      mnar_data_sub_wgt <- mnar_data_sub$wgt

      # imputation (no prs)
      these_vars    <- c(outcome, exposure, covs)
      mcar_imp_data <- mice(mcar_data_sub |> dplyr::select(all_of(these_vars)), printFlag = FALSE)
      mar_imp_data  <- mice(mar_data_sub |> dplyr::select(all_of(these_vars)), printFlag = FALSE)
      mnar_imp_data <- mice(mnar_data_sub |> dplyr::select(all_of(these_vars)), printFlag = FALSE)

      mcar_dsn <- svydesign(
        ids     = ~1,
        weights = ~wgt,
        data    = mcar_data_sub
      )
      mar_dsn <- svydesign(
        ids     = ~1,
        weights = ~wgt,
        data    = mar_data_sub
      )
      mnar_dsn <- svydesign(
        ids     = ~1,
        weights = ~wgt,
        data    = mnar_data_sub
      )

      # imputation w/ prs
      these_vars_prs    <- c(these_vars, "bmi_prs", "glucose_prs")
      mcar_imp_data_prs <- mice(mcar_data_sub |> dplyr::select(all_of(these_vars_prs)), printFlag = FALSE)
      mar_imp_data_prs  <- mice(mar_data_sub |> dplyr::select(all_of(these_vars_prs)), printFlag = FALSE)
      mnar_imp_data_prs <- mice(mnar_data_sub |> dplyr::select(all_of(these_vars_prs)), printFlag = FALSE)

      rbindlist(
        list(
          # truth
          .glm(as.formula(f), data = data, var = exposure)[, `:=`(method = "Full", adj = "Unadjusted", miss_data = NA_character_, weight = NA_character_)],
          .glm(as.formula(f_cov), data = data, var = exposure)[, `:=`(method = "Full", adj = "Adjusted", miss_data = NA_character_, weight = NA_character_)],

          # mcar
          ## unadjusted
          .glm(as.formula(f), data = mcar_data_sub, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MCAR", weight = "Unweighted")],
          .glm_pool(as.formula(f), mcar_imp_data, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MCAR", weight = "Unweighted")],
          .glm_pool(as.formula(f), mcar_imp_data_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MCAR", weight = "Unweighted")],
          .svyglm(as.formula(f), design = mcar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MCAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f), imp_object = mcar_imp_data, weights = mcar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MCAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f), imp_object = mcar_imp_data_prs, weights = mcar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MCAR", weight = "Weighted")],
          ## adjusted
          .glm(as.formula(f_cov), data = mcar_data_sub, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MCAR", weight = "Unweighted")],
          .glm_pool(as.formula(f_cov), mcar_imp_data, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MCAR", weight = "Unweighted")],
          .glm_pool(as.formula(f_cov), mcar_imp_data_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MCAR", weight = "Unweighted")],
          .svyglm(as.formula(f_cov), design = mcar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MCAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f_cov), imp_object = mcar_imp_data, weights = mcar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MCAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f_cov), imp_object = mcar_imp_data_prs, weights = mcar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MCAR", weight = "Weighted")],

          # mar
          ## unadjusted
          .glm(as.formula(f), data = mar_data_sub, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MAR", weight = "Unweighted")],
          .glm_pool(as.formula(f), mar_imp_data, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MAR", weight = "Unweighted")],
          .glm_pool(as.formula(f), mar_imp_data_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MAR", weight = "Unweighted")],
          .svyglm(as.formula(f), design = mar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f), imp_object = mar_imp_data, weights = mar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f), imp_object = mar_imp_data_prs, weights = mar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MAR", weight = "Weighted")],
          ## adjusted
          .glm(as.formula(f_cov), data = mar_data_sub, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MAR", weight = "Unweighted")],
          .glm_pool(as.formula(f_cov), mar_imp_data, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MAR", weight = "Unweighted")],
          .glm_pool(as.formula(f_cov), mar_imp_data_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MAR", weight = "Unweighted")],
          .svyglm(as.formula(f_cov), design = mar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f_cov), imp_object = mar_imp_data, weights = mar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f_cov), imp_object = mar_imp_data_prs, weights = mar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MAR", weight = "Weighted")],

          #mnar
          ## unadjusted
          .glm(as.formula(f), data = mnar_data_sub, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MNAR", weight = "Unweighted")],
          .glm_pool(as.formula(f), mnar_imp_data, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MNAR", weight = "Unweighted")],
          .glm_pool(as.formula(f), mnar_imp_data_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MNAR", weight = "Unweighted")],
          .svyglm(as.formula(f), design = mnar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Unadjusted", miss_data = "MNAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f), imp_object = mnar_imp_data, weights = mnar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed", adj = "Unadjusted", miss_data = "MNAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f), imp_object = mnar_imp_data_prs, weights = mnar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Unadjusted", miss_data = "MNAR", weight = "Weighted")],
          ## adjusted
          .glm(as.formula(f_cov), data = mnar_data_sub, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MNAR", weight = "Unweighted")],
          .glm_pool(as.formula(f_cov), mnar_imp_data, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MNAR", weight = "Unweighted")],
          .glm_pool(as.formula(f_cov), mnar_imp_data_prs, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MNAR", weight = "Unweighted")],
          .svyglm(as.formula(f_cov), design = mnar_dsn, var = exposure)[, `:=`(method = "Complete case", adj = "Adjusted", miss_data = "MNAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f_cov), imp_object = mnar_imp_data, weights = mnar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed", adj = "Adjusted", miss_data = "MNAR", weight = "Weighted")],
          .svyglm_pool(as.formula(f_cov), imp_object = mnar_imp_data_prs, weights = mnar_data_sub_wgt, var = exposure)[, `:=`(method = "Imputed w/ PRS", adj = "Adjusted", miss_data = "MNAR", weight = "Weighted")]
        ),
        fill = TRUE, use.names = TRUE
      )[, `:=` (seed = i, size = sample_size)]
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )

  write_csv(
    x    = tmp,
    file = paste0("data/public/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_res.csv")
  )

  #####
  truth_una <- tmp[method == "Full" & adj == "Unadjusted", ]
  truth_adj <- tmp[method == "Full" & adj == "Adjusted", ]

  mcar_cc_unw_una   <- tmp[method == "Complete case" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mcar_imp_unw_una  <- tmp[method == "Imputed" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mcar_impp_unw_una <- tmp[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mcar_cc_w_una     <- tmp[method == "Complete case" & miss_data == "MCAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mcar_imp_w_una    <- tmp[method == "Imputed" & miss_data == "MCAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mcar_impp_w_una   <- tmp[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Weighted" & adj == "Unadjusted", ]

  mcar_cc_unw_adj   <- tmp[method == "Complete case" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mcar_imp_unw_adj  <- tmp[method == "Imputed" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mcar_impp_unw_adj <- tmp[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mcar_cc_w_adj     <- tmp[method == "Complete case" & miss_data == "MCAR" & weight == "Weighted" & adj == "Adjusted", ]
  mcar_imp_w_adj    <- tmp[method == "Imputed" & miss_data == "MCAR" & weight == "Weighted" & adj == "Adjusted", ]
  mcar_impp_w_adj   <- tmp[method == "Imputed w/ PRS" & miss_data == "MCAR" & weight == "Weighted" & adj == "Adjusted", ]

  mar_cc_unw_una   <- tmp[method == "Complete case" & miss_data == "MAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mar_imp_unw_una  <- tmp[method == "Imputed" & miss_data == "MAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mar_impp_unw_una <- tmp[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mar_cc_w_una     <- tmp[method == "Complete case" & miss_data == "MAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mar_imp_w_una    <- tmp[method == "Imputed" & miss_data == "MAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mar_impp_w_una   <- tmp[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Weighted" & adj == "Unadjusted", ]

  mar_cc_unw_adj   <- tmp[method == "Complete case" & miss_data == "MAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mar_imp_unw_adj  <- tmp[method == "Imputed" & miss_data == "MAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mar_impp_unw_adj <- tmp[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mar_cc_w_adj     <- tmp[method == "Complete case" & miss_data == "MAR" & weight == "Weighted" & adj == "Adjusted", ]
  mar_imp_w_adj    <- tmp[method == "Imputed" & miss_data == "MAR" & weight == "Weighted" & adj == "Adjusted", ]
  mar_impp_w_adj   <- tmp[method == "Imputed w/ PRS" & miss_data == "MAR" & weight == "Weighted" & adj == "Adjusted", ]

  mnar_cc_unw_una   <- tmp[method == "Complete case" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mnar_imp_unw_una  <- tmp[method == "Imputed" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mnar_impp_unw_una <- tmp[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Unadjusted", ]
  mnar_cc_w_una     <- tmp[method == "Complete case" & miss_data == "MNAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mnar_imp_w_una    <- tmp[method == "Imputed" & miss_data == "MNAR" & weight == "Weighted" & adj == "Unadjusted", ]
  mnar_impp_w_una   <- tmp[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Weighted" & adj == "Unadjusted", ]

  mnar_cc_unw_adj   <- tmp[method == "Complete case" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mnar_imp_unw_adj  <- tmp[method == "Imputed" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mnar_impp_unw_adj <- tmp[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Unweighted" & adj == "Adjusted", ]
  mnar_cc_w_adj     <- tmp[method == "Complete case" & miss_data == "MNAR" & weight == "Weighted" & adj == "Adjusted", ]
  mnar_imp_w_adj    <- tmp[method == "Imputed" & miss_data == "MNAR" & weight == "Weighted" & adj == "Adjusted", ]
  mnar_impp_w_adj   <- tmp[method == "Imputed w/ PRS" & miss_data == "MNAR" & weight == "Weighted" & adj == "Adjusted", ]

  simul_diag <- function(cleaned, truth) {
    raw_bias <- mean(cleaned$est - truth)
    per_bias <- mean(100 * abs(raw_bias / truth))
    cov_rate <- mean(cleaned$lo < truth & cleaned$hi > truth)
    avg_wide <- mean(cleaned$hi - cleaned$lo)
    rmse     <- sqrt(mean((cleaned$est - truth)^2))
    data.table(
      adj       = cleaned[, unique(adj)],
      miss_data = cleaned[, unique(miss_data)],
      weight    = cleaned[, unique(weight)],
      size      = cleaned[, unique(size)],
      method    = cleaned[, unique(method)],
      raw_bias  = raw_bias,
      per_bias  = per_bias,
      cov_rate  = cov_rate,
      avg_wide  = avg_wide,
      rmse      = rmse
    )
  }

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

  write_csv(
    x = diag,
    file = paste0("data/public/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_diag.csv")
  )
}


# diag |>
#   dplyr::mutate(
#     adj = factor(adj, levels = c("Unadjusted", "Adjusted")),
#     miss_dat = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
#     # analysis = factor(analysis, levels = c("Sample (CC)", "Sample (IMP)", "Sample (IMP w/ PRS)", "Weighted sample (CC)", "Weighted sample (IMP w/ PRS)"))
#     method = factor(method, levels = c("Complete case", "Imputed", "Imputed w/ PRS")),
#     weight = factor(weight, levels = c("Unweighted", "Weighted"))
#   ) |>
#   dplyr::select(size, adj, method, miss_dat, weight, cov_rate) |>
#   ggplot(aes(x = size, y = cov_rate, shape = method, color = miss_dat, linetype = method)) +
#   geom_hline(yintercept = 0.95, linetype = "dashed") +
#   # geom_line(linewidth = 1) +
#   geom_point(size = 3) +
#   scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(
#     title = stringr::str_wrap("Coverage rate of 95% CI for BMI coefficient for glucose by missing data mechanism and method", 75),
#     x = "Sample size",
#     y = "Coverage rate (%)",
#     caption = "Dashed line represents 95% coverage rate"
#   ) +
#   facet_wrap(~adj+weight) +
#   scale_color_ms() +
#   theme_ms() +
#   theme(
#     legend.title = element_blank()
#   )

# ####
# map_generate_data <- function(
#   iterations = 1000,
#   size = 50000,
#   mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
#   mat = NULL
# ) {
#   purrr::map(
#     seq_len(iterations),
#     \(x) {
#       set.seed(x)
#       generate_data(n = size, mu = mu, mat = mat)
#     },
#     .progress = "generating data")
# }

# sample_data <- function(
#   data_list,
#   size,
#   prob = NULL,
#   reps = 1,
#   replace = TRUE,
#   boot = FALSE
# )

# sample_sizes <- c(1000, 2500, 5000)
# iterations   <- 1000
# outcome      <- "glucose"
# exposure     <- "bmi"
# covs         <- c("age", "female", "nhw", "smoke")
# f            <- paste0(outcome, " ~ ", exposure)
# f_cov        <- paste0(f, " + ", paste(covs, collapse = " + "))

# source("fn/resample_helpers.R")
# source("fn/regression_helpers.R")

# formals(map_dfr)$.progress <- formals(map)$.progress <- TRUE
# formals(mice)$printFlag <- FALSE


#   # generate dataset
#     datas <- map_generate_data(iterations = 1, size = 1000, mat = mat)
#   # sample from dataset


# for (sample_size in sample_sizes) {
#   print(sample_size)

#   samples <- map_generate_data(iterations = iterations, size = sample_size, mat = mat)
  
#   # truth
#   true_un_betas <- map_dfr(samples, \(x) .glm(f, data = x, var = exposure))
#   true_adj_betas <- map_dfr(samples, \(x) .glm(f_cov, data = x, var = exposure))
  
#   mcar_dat <- map(
#     samples,
#     \(x) {
#       mcar_ind <- rbinom(nrow(x), 1, 0.25)
#       x[mcar_ind == 1, exposure] <- NA
#       x
#     }
#   )
#   mar_dat <- map(
#     samples,
#     \(x) {
#       cbind(x[, .(id)], make_missing(
#         data       = x,
#         target_var = exposure,
#         mech       = "MAR",
#         vars       = c(outcome, covs),
#         scale_cont = FALSE
#       ))
#     }
#   )
#   mnar_dat <- map(
#     samples,
#     \(x) {
#       cbind(x[, .(id)], make_missing(
#         data       = x,
#         target_var = exposure,
#         mech       = "MNAR",
#         vars       = c(outcome, covs),
#         intercept  = -2.75,
#         scale_cont = FALSE
#       ))
#     }
#   )
  
#   # imputation
#   cli_progress_step("Imputing data")
#   these_vars <- c(outcome, exposure, covs)
#   mcar_imp <- map(mcar_dat, \(x) mice(data = x[, ..these_vars]))
#   mar_imp  <- map(mar_dat, \(x) mice(data = x[, ..these_vars]))
#   mnar_imp <- map(mnar_dat, \(x) mice(data = x[, ..these_vars]))
  
#   # imputation with prs
#   cli_progress_step("Imputing data with PRS")
#   these_vars_prs <- c(these_vars, "bmi_prs", "glucose_prs")
#   mcar_imp_prs <- map(mcar_dat, \(x) mice(x[, ..these_vars_prs]))
#   mar_imp_prs <- map(mar_dat, \(x) mice(x[, ..these_vars_prs]))
#   mnar_imp_prs <- map(mnar_dat, \(x) mice(x[, ..these_vars_prs]))
  
#   # List of data sets
#   data_sets <- list(
#     mcar_dat     = mcar_dat,
#     mar_dat      = mar_dat,
#     mnar_dat     = mnar_dat,
#     mcar_imp     = mcar_imp,
#     mar_imp      = mar_imp,
#     mnar_imp     = mnar_imp,
#     mcar_imp_prs = mcar_imp_prs,
#     mar_imp_prs  = mar_imp_prs,
#     mnar_imp_prs = mnar_imp_prs
#   )
  
#   # List to store results
#   results <- list()
  
#   # Loop through data sets and apply .glm_pool for both models
#   for (data_name in names(data_sets)) {
#     cli::cli_progress_step(paste0("fitting models: ", data_name))
#     un_model_name <- paste0(data_name, "_un_betas")
#     adj_model_name <- paste0(data_name, "_adj_betas")
#     if ("mids" %in% class(data_sets[[data_name]][[1]])) {
#       results[[un_model_name]]  <- map_glm_pool(data_sets[[data_name]], f, exposure)
#       results[[adj_model_name]] <- map_glm_pool(data_sets[[data_name]], f_cov, exposure)
#     } else {
#       results[[un_model_name]]  <- map_glm(data_sets[[data_name]], f, exposure)
#       results[[adj_model_name]] <- map_glm(data_sets[[data_name]], f_cov, exposure)
#     }
#   }
  
#   # summarize
#   means <- tribble(
#     ~adj, ~est, ~method, ~miss_dat,
#     "Unadjusted", mean(true_un_betas$est), "True", NA,
#     "Unadjusted", mean(results[["mcar_dat_un_betas"]]$est), "Complete Case", "MCAR",
#     "Unadjusted", mean(results[["mar_dat_un_betas"]]$est), "Complete Case", "MAR",
#     "Unadjusted", mean(results[["mnar_dat_un_betas"]]$est), "Complete Case", "MNAR",
#     "Unadjusted", mean(results[["mcar_imp_un_betas"]]$est), "Multiple Imputation", "MCAR",
#     "Unadjusted", mean(results[["mar_imp_un_betas"]]$est), "Multiple Imputation", "MAR",
#     "Unadjusted", mean(results[["mnar_imp_un_betas"]]$est), "Multiple Imputation", "MNAR",
#     "Unadjusted", mean(results[["mcar_imp_prs_un_betas"]]$est), "Multiple Imputation w/ PRS", "MCAR",
#     "Unadjusted", mean(results[["mar_imp_prs_un_betas"]]$est), "Multiple Imputation w/ PRS", "MAR",
#     "Unadjusted", mean(results[["mnar_imp_prs_un_betas"]]$est), "Multiple Imputation w/ PRS", "MNAR",
#     "Adjusted", mean(true_adj_betas$est), "True", NA,
#     "Adjusted", mean(results[["mcar_dat_adj_betas"]]$est), "Complete Case", "MCAR",
#     "Adjusted", mean(results[["mar_dat_adj_betas"]]$est), "Complete Case", "MAR",
#     "Adjusted", mean(results[["mnar_dat_adj_betas"]]$est), "Complete Case", "MNAR",
#     "Adjusted", mean(results[["mcar_imp_adj_betas"]]$est), "Multiple Imputation", "MCAR",
#     "Adjusted", mean(results[["mar_imp_adj_betas"]]$est), "Multiple Imputation", "MAR",
#     "Adjusted", mean(results[["mnar_imp_adj_betas"]]$est), "Multiple Imputation", "MNAR",
#     "Adjusted", mean(results[["mcar_imp_prs_adj_betas"]]$est), "Multiple Imputation w/ PRS", "MCAR",
#     "Adjusted", mean(results[["mar_imp_prs_adj_betas"]]$est), "Multiple Imputation w/ PRS", "MAR",
#     "Adjusted", mean(results[["mnar_imp_prs_adj_betas"]]$est), "Multiple Imputation w/ PRS", "MNAR"
#   )
  
#   diag_tab <- bind_rows(
#     cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MCAR"), simul_diag(results[["mcar_dat_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MAR"), simul_diag(results[["mar_dat_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MNAR"), simul_diag(results[["mnar_dat_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MAR"), simul_diag(results[["mar_imp_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_prs_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MAR"), simul_diag(results[["mar_imp_prs_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_prs_un_betas"]], mean(true_un_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MCAR"), simul_diag(results[["mcar_dat_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MAR"), simul_diag(results[["mar_dat_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MNAR"), simul_diag(results[["mnar_dat_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MAR"), simul_diag(results[["mar_imp_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_prs_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MAR"), simul_diag(results[["mar_imp_prs_adj_betas"]], mean(true_adj_betas$est))),
#     cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_prs_adj_betas"]], mean(true_adj_betas$est)))
#   )
  
#   write_delim(
#     means,
#     glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_means.csv"),
#     delim = ","
#   )
#   write_delim(
#     diag_tab,
#     glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_diag.csv"),
#     delim = ","
#   )
#   qsave(
#     results,
#     glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_results.qs")
#   )
# }
