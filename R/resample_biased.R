# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr, glue, cli, glue
)

plan(multicore, workers = 10)

sample_sizes <- c(1000, 2500, 5000)
iterations   <- 100
outcome      <- "glucose"
exposure     <- "bmi"
covs         <- c("age", "female", "nhw", "smoke")
f            <- paste0(outcome, " ~ ", exposure)
f_cov        <- paste0(f, " + ", paste(covs, collapse = " + "))

# data -------------------------------------------------------------------------
analytic_sub <- qread(glue("data/private/analytic_sub.qs"))
pim          <- qread(glue("data/private/mgi/20220822/MGI_PIM0X_20220822.qs"))
pim          <- pim[, .(id, t2d = EM_202.2)]
prs          <- qread("data/private/ex_prs_processed.qs")
prs          <- prs[, .(id, bmi_prs = bmi, glucose_prs = glucose)]

# merge pim vars into analytic_sub
merged <- analytic_sub |> left_join(pim, by = "id")

# keep limited data
d <- merged |>
    mutate(nhw = as.numeric(race_eth == "NHW")) |>
    select(
        id,
        age, female,
        bmi = bmi_med, smoke = smoke_ever,
        glucose = glucose_med, nhw
    )
d_id <- d |> select(id)
# d <- d |> select(-id)

data <- d |> select(id, age, female, glucose, bmi, smoke, nhw)

scaled <- data |>
    mutate(
        age     = scale(age)[, 1],
        glucose = scale(glucose)[, 1],
        bmi     = scale(bmi)[, 1]
    )

scaled <- inner_join(
    scaled,
    prs,
    by = "id"
)

cor(scaled[, !c("id")])

# functions --------------------------------------------------------------------
source("fn/resample_helpers.R")
source("fn/regression_helpers.R")

formals(map_dfr)$.progress <- formals(map)$.progress <- TRUE
formals(mice)$printFlag <- FALSE

# create weights
weight_score <- scaled$age + scaled$glucose + scaled$bmi
weight_prob <- plogis(weight_score)

scaled$wgt_prob <- weight_prob
scaled$wgt      <- 1 / weight_prob
scaled$wgt <- scaled$wgt * (nrow(scaled) / sum(scaled$wgt))

true_un_beta <- .glm(f, data = scaled, var = exposure)$est
true_adj_beta <- .glm(f_cov, data = scaled, var = exposure)$est

# analysis ---------------------------------------------------------------------
for (sample_size in sample_sizes) {
    print(sample_size)
    # create samples with weights
    samples <- re_sample(scaled, sample_size, iterations, prob = scaled$wgt_prob)
    samples_wgts <- map(samples, \(x) {
        x[, wgt] * (nrow(x) / sum(x[, wgt]))
    })   
    # truth
    # true_un_betas <- map_dfr(samples, \(x) .glm(f, data = x, var = exposure))
    # true_adj_betas <- map_dfr(samples, \(x) .glm(f_cov, data = x, var = exposure))

    mcar_dat <- map(
        samples,
        \(x) {
            mcar_ind <- rbinom(nrow(x), 1, 0.25)
            x[mcar_ind == 1, exposure] <- NA
            x
        }
    )
    mar_dat <- map(
        samples,
        \(x) {
            cbind(x[, .(id)], make_missing(
                data       = x[, !c("id")],
                target_var = exposure,
                vars       = covs,
                mech       = "MAR",
                scale_cont = FALSE
            ))
        }
    )
    mnar_dat <- map(
        samples,
        \(x) {
            cbind(x[, .(id)], make_missing(
                data       = x[, !c("id")],
                target_var = exposure,
                mech       = "MNAR",
                intercept  = -2.75,
                scale_cont = FALSE
            ))
        }
    )

    # imputation
    cli_progress_step("Imputing data")
    these_vars <- c(outcome, exposure, covs)
    mcar_imp <- map(mcar_dat, \(x) mice(data = x[, ..these_vars]))
    mar_imp  <- map(mar_dat, \(x) mice(data = x[, ..these_vars]))
    mnar_imp <- map(mnar_dat, \(x) mice(data = x[, ..these_vars]))

    # imputation with prs
    cli_progress_step("Imputing data with PRS")
    these_vars_prs <- c(outcome, exposure, covs, "bmi_prs", "glucose_prs")
    mcar_imp_prs <- map(mcar_dat, \(x) mice(merge(x, prs, "id")[, ..these_vars_prs]))
    mar_imp_prs  <- map(mar_dat, \(x) mice(merge(x, prs, "id")[, ..these_vars_prs]))
    mnar_imp_prs <- map(mnar_dat, \(x) mice(merge(x, prs, "id")[, ..these_vars_prs]))

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
        wun_model_name <- paste0(data_name, "_wun_betas")
        adj_model_name <- paste0(data_name, "_adj_betas")
        wadj_model_name <- paste0(data_name, "_wadj_betas")
        if ("mids" %in% class(data_sets[[data_name]][[1]])) {
            results[[un_model_name]]  <- map_glm_pool(data_sets[[data_name]], f, exposure)
            results[[adj_model_name]] <- map_glm_pool(data_sets[[data_name]], f_cov, exposure)
            results[[wun_model_name]] <- map_svyglm_pool(data_sets[[data_name]], formula = f, weights_list = samples_wgts, var = exposure)
            results[[wadj_model_name]] <- map_svyglm_pool(data_sets[[data_name]], formula = f_cov, weights_list = samples_wgts, var = exposure)
        } else {
            results[[un_model_name]]  <- map_glm(data_sets[[data_name]], f, exposure)
            results[[adj_model_name]] <- map_glm(data_sets[[data_name]], f_cov, exposure)
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
        cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MCAR"), simul_diag(results[["mcar_dat_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MAR"), simul_diag(results[["mar_dat_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Complete case", miss_dat = "MNAR"), simul_diag(results[["mnar_dat_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MAR"), simul_diag(results[["mar_imp_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_prs_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MAR"), simul_diag(results[["mar_imp_prs_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Unadjusted", method = "Imputed w/ PRS", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_prs_un_betas"]], mean(true_un_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MCAR"), simul_diag(results[["mcar_dat_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MAR"), simul_diag(results[["mar_dat_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Complete case", miss_dat = "MNAR"), simul_diag(results[["mnar_dat_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MAR"), simul_diag(results[["mar_imp_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MCAR"), simul_diag(results[["mcar_imp_prs_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MAR"), simul_diag(results[["mar_imp_prs_adj_betas"]], mean(true_adj_betas$est))),
        cbind(tibble(adj = "Adjusted", method = "Imputed w/ PRS", miss_dat = "MNAR"), simul_diag(results[["mnar_imp_prs_adj_betas"]], mean(true_adj_betas$est)))
    )

    write_delim(
        means,
        glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_means.csv"),
        delim = ","
    )
    write_delim(
        diag_tab,
        glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_diag.csv"),
        delim = ","
    )
    qsave(
        results,
        glue("results/{outcome}_{exposure}_n{sample_size}_i{iterations}_results.qs")
    )
}
