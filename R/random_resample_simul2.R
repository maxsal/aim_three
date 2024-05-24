# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr, glue, cli, glue
)

if (availableCores() <= 8) {
    plan(multisession, workers = 6)
} else {
    wkrs <- min(c(availableCores() / 2, 24))
    plan(multicore, workers = wkrs)
}

outcome      <- "glucose"
exposure     <- "bmi"
covs         <- c("age", "female", "nhw", "smoke")
sample_sizes <- c(1000, 2500, 5000, 10000)
iterations   <- 1000
f            <- paste0(outcome, " ~ ", exposure)
f_cov        <- paste0(f, " + ", paste0(covs, collapse = " + "))
ms <- 5 # number of imputations

# data -------------------------------------------------------------------------
means <- c(
    glucose = -1.96825433778349e-16, bmi = -7.91723213769746e-17,
    age = 1.20111356353818e-16, female = 1.42619696156583e-16, nhw = -5.03312408648097e-18,
    smoke = -3.2111110600884e-17, bmi_prs = -1.15519592499013e-14,
    glucose_prs = -7.07320752607373e-15
)

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
    n = 100000,
    mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
    bin_vars = c("female", "smoke", "nhw"),
    mat = NULL) {
    data <- MASS::mvrnorm(
        n = n,
        mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
        Sigma = mat
    ) |> data.table::as.data.table()
    # data <- mvnfast::rmvn(
    #     n = n,
    #     mu = mu,
    #     sigma = mat
    # ) |> data.table::as.data.table()
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

# functions --------------------------------------------------------------------
source("fn/resample_helpers.R")
source("fn/regression_helpers.R")

formals(map_dfr)$.progress   <- formals(future_map)$.progress <- formals(future_map_dfr)$.progress <- formals(map)$.progress <- TRUE
formals(future_map)$.options <- formals(future_map_dfr)$.options <- furrr_options(seed = TRUE)
formals(mice)$printFlag      <- FALSE

# analysis ---------------------------------------------------------------------
for (sample_size in sample_sizes) {
    print(sample_size)

    cli_progress_step("generating data...")
    datas   <- map_generate_data(mat = mat, mu = means, iterations = iterations)
    samples <- map(datas, \(x) re_sample(x = x, size = sample_size, reps = 1)[[1]])

    # truth
    true_un_betas  <- future_map_dfr(datas, \(x) .glm(f, data = x, var = exposure), .progress = TRUE, .options = furrr_options(seed = TRUE))
    true_adj_betas <- future_map_dfr(datas, \(x) .glm(f_cov, data = x, var = exposure), .progress = TRUE, .options = furrr_options(seed = TRUE))

    mcar_dat <- map(
        samples,
        \(x) {
            mcar_ex_ind  <- rbinom(nrow(x), 1, 0.25)
            mcar_out_ind <- rbinom(nrow(x), 1, 0.25)
            x[mcar_ex_ind == 1, exposure] <- NA
            x[mcar_out_ind == 1, outcome] <- NA
            x
        }
    )

    mar_dat <- map(
        samples,
        \(x) {
            make_ex_out_mar5(
                data = x,
                exposure = exposure,
                outcome = outcome,
                covs = covs
            )
        }
    )
    
    mnar_dat <- map(
        samples,
        \(x) {
            make_ex_out_mnar5(
                data = x,
                exposure = exposure,
                outcome = outcome,
                covs = covs
            )
        }
    )

    # imputation
    cli_progress_step("Imputing data")
    these_vars <- c(outcome, exposure, covs)
    
    mcar_imp <- future_map(mcar_dat, \(x) mice(data = x |> select(all_of(these_vars)), m = ms), .progress = TRUE, .options = furrr_options(seed = TRUE))
    mar_imp <- future_map(mar_dat, \(x) mice(data = x |> select(all_of(these_vars)), m = ms), .progress = TRUE, .options = furrr_options(seed = TRUE))
    mnar_imp <- future_map(mnar_dat, \(x) mice(data = x |> select(all_of(these_vars)), m = ms), .progress = TRUE, .options = furrr_options(seed = TRUE))

    # imputation with prs
    cli_progress_step("Imputing data with PRS")
    these_vars_prs <- c(these_vars, "bmi_prs", "glucose_prs")
    mcar_imp_prs <- future_map(mcar_dat, \(x) mice(x |> select(all_of(these_vars)), m = ms), .progress = TRUE, .options = furrr_options(seed = TRUE))
    mar_imp_prs <- future_map(mar_dat, \(x) mice(x |> select(all_of(these_vars)), m = ms), .progress = TRUE, .options = furrr_options(seed = TRUE))
    mnar_imp_prs <- future_map(mnar_dat, \(x) mice(x |> select(all_of(these_vars)), m = ms), .progress = TRUE, .options = furrr_options(seed = TRUE))

    # List of data sets
    data_sets <- list(
        mcar_dat     = mcar_dat,
        mar_dat      = mar_dat,
        mnar_dat     = mnar_dat,
        mcar_imp     = mcar_imp,
        mar_imp      = mar_imp,
        mnar_imp     = mnar_imp,
        mcar_imp_prs = mcar_imp_prs,
        mar_imp_prs  = mar_imp_prs,
        mnar_imp_prs = mnar_imp_prs
    )

    # List to store results
    results <- list()

    # Loop through data sets and apply .glm_pool for both models
    for (data_name in names(data_sets)) {
        cli::cli_progress_step(paste0("fitting models: ", data_name))
        un_model_name <- paste0(data_name, "_un_betas")
        adj_model_name <- paste0(data_name, "_adj_betas")
        if ("mids" %in% class(data_sets[[data_name]][[1]])) {
            results[[un_model_name]]  <- future_map_glm_pool(data_sets[[data_name]], f, exposure)
            results[[adj_model_name]] <- future_map_glm_pool(data_sets[[data_name]], f_cov, exposure)
        } else {
            results[[un_model_name]]  <- future_map_glm(data_sets[[data_name]], f, exposure)
            results[[adj_model_name]] <- future_map_glm(data_sets[[data_name]], f_cov, exposure)
        }
    }

    # summarize
    means <- tribble(
        ~adj, ~est, ~method, ~miss_dat,
        "Unadjusted", mean(true_un_betas$est), "True", NA,
        "Unadjusted", mean(results[["mcar_dat_un_betas"]]$est), "Complete Case", "MCAR",
        "Unadjusted", mean(results[["mar_dat_un_betas"]]$est), "Complete Case", "MAR",
        "Unadjusted", mean(results[["mnar_dat_un_betas"]]$est), "Complete Case", "MNAR",
        "Unadjusted", mean(results[["mcar_imp_un_betas"]]$est), "Multiple Imputation", "MCAR",
        "Unadjusted", mean(results[["mar_imp_un_betas"]]$est), "Multiple Imputation", "MAR",
        "Unadjusted", mean(results[["mnar_imp_un_betas"]]$est), "Multiple Imputation", "MNAR",
        "Unadjusted", mean(results[["mcar_imp_prs_un_betas"]]$est), "Multiple Imputation w/ PRS", "MCAR",
        "Unadjusted", mean(results[["mar_imp_prs_un_betas"]]$est), "Multiple Imputation w/ PRS", "MAR",
        "Unadjusted", mean(results[["mnar_imp_prs_un_betas"]]$est), "Multiple Imputation w/ PRS", "MNAR",
        "Adjusted", mean(true_adj_betas$est), "True", NA,
        "Adjusted", mean(results[["mcar_dat_adj_betas"]]$est), "Complete Case", "MCAR",
        "Adjusted", mean(results[["mar_dat_adj_betas"]]$est), "Complete Case", "MAR",
        "Adjusted", mean(results[["mnar_dat_adj_betas"]]$est), "Complete Case", "MNAR",
        "Adjusted", mean(results[["mcar_imp_adj_betas"]]$est), "Multiple Imputation", "MCAR",
        "Adjusted", mean(results[["mar_imp_adj_betas"]]$est), "Multiple Imputation", "MAR",
        "Adjusted", mean(results[["mnar_imp_adj_betas"]]$est), "Multiple Imputation", "MNAR",
        "Adjusted", mean(results[["mcar_imp_prs_adj_betas"]]$est), "Multiple Imputation w/ PRS", "MCAR",
        "Adjusted", mean(results[["mar_imp_prs_adj_betas"]]$est), "Multiple Imputation w/ PRS", "MAR",
        "Adjusted", mean(results[["mnar_imp_prs_adj_betas"]]$est), "Multiple Imputation w/ PRS", "MNAR"
    )

    diag_tab <- bind_rows(
        cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MCAR"), simul_diag2(results[["mcar_dat_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MAR"), simul_diag2(results[["mar_dat_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MNAR"), simul_diag2(results[["mnar_dat_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MCAR"), simul_diag2(results[["mcar_imp_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MAR"), simul_diag2(results[["mar_imp_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MNAR"), simul_diag2(results[["mnar_imp_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MCAR"), simul_diag2(results[["mcar_imp_prs_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MAR"), simul_diag2(results[["mar_imp_prs_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MNAR"), simul_diag2(results[["mnar_imp_prs_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MCAR"), simul_diag2(results[["mcar_dat_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MAR"), simul_diag2(results[["mar_dat_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MNAR"), simul_diag2(results[["mnar_dat_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MCAR"), simul_diag2(results[["mcar_imp_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MAR"), simul_diag2(results[["mar_imp_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MNAR"), simul_diag2(results[["mnar_imp_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MCAR"), simul_diag2(results[["mcar_imp_prs_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MAR"), simul_diag2(results[["mar_imp_prs_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MNAR"), simul_diag2(results[["mnar_imp_prs_adj_betas"]], mean(true_adj_betas$est)))
    )

    write_delim(
        means,
        glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_m{ms}_exoutmiss_means_random.csv"),
        delim = ","
    )
    write_delim(
        diag_tab,
        glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_m{ms}_exoutmiss_diag_random.csv"),
        delim = ","
    )
    qsave(
        results,
        glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_m{ms}_exoutmiss_results_random.qs")
    )
}
