# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr, glue, cli, glue
)

if (availableCores() <= 8) {
    plan(multisession, workers = 6)
} else if (availableCores() >= 64) {
    plan(multicore, workers = 18)
}


# sample_sizes <- c(1000, 2500, 5000, 10000)
# iterations   <- 1000
# outcome      <- "glucose"
# exposure     <- "bmi"
# covs         <- c("age", "female", "nhw", "smoke")
# f            <- paste0(outcome, " ~ ", exposure)
# f_cov        <- paste0(f, " + ", paste(covs, collapse = " + "))

outcome      <- "glucose"
exposure     <- "bmi"
covs         <- c("age", "female", "nhw", "smoke")
sample_sizes <- c(1000, 2500, 5000, 10000)
iterations   <- 1000
f            <- paste0(outcome, " ~ ", exposure)
f_cov        <- paste0(f, " + ", paste0(covs, collapse = " + "))

# data -------------------------------------------------------------------------
cor_structure <- structure(
    c(
        1, -0.138472771069173, 0.262539473913268, 0.0468405451550121,
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
        1
    ),
    dim = c(8L, 8L),
    dimnames = list(
        c(
            "age", "female", "glucose", "bmi", "smoke", "nhw",
            "bmi_prs", "glucose_prs"
        ),
        c(
            "age", "female", "glucose", "bmi", "smoke", "nhw",
            "bmi_prs", "glucose_prs"
        )
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
    n = 100000,
    mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
    bin_vars = c("female", "smoke", "nhw"),
    mat = NULL) {
    data <- MASS::mvrnorm(
        n = n,
        mu = c(0, 0.5359, 0, 0, 0.4811, 0.8864, 0, 0),
        Sigma = mat
    ) |> data.table::as.data.table()
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
    purrr::map(
        seq_len(iterations),
        \(x) {
            set.seed(x)
            generate_data(n = size, mu = mu, mat = mat)
        },
        .progress = "generating data"
    )
}

# functions --------------------------------------------------------------------
source("fn/resample_helpers.R")
source("fn/regression_helpers.R")

formals(map_dfr)$.progress <- formals(map)$.progress <- TRUE
formals(mice)$printFlag <- FALSE

scaled <- generate_data(mat = mat)

# analysis ---------------------------------------------------------------------
for (sample_size in sample_sizes) {
    print(sample_size)

    samples <- re_sample(scaled, sample_size, iterations)
        
    # truth
    true_un_betas <- map_dfr(samples, \(x) .glm(f, data = x, var = exposure))
    true_adj_betas <- map_dfr(samples, \(x) .glm(f_cov, data = x, var = exposure))

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
    mcar_imp <- map(mcar_dat, \(x) mice(data = x[, !c("id")]))
    mar_imp  <- map(mar_dat, \(x) mice(data = x[, !c("id")]))
    mnar_imp <- map(mnar_dat, \(x) mice(data = x[, !c("id")]))

    # imputation with prs
    cli_progress_step("Imputing data with PRS")
    mcar_imp_prs <- map(mcar_dat, \(x) mice(merge(x, prs, "id")[, !c("id")]))
    mar_imp_prs  <- map(mar_dat, \(x) mice(merge(x, prs, "id")[, !c("id")]))
    mnar_imp_prs <- map(mnar_dat, \(x) mice(merge(x, prs, "id")[, !c("id")]))

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
            results[[un_model_name]]  <- map_glm_pool(data_sets[[data_name]], f, exposure)
            results[[adj_model_name]] <- map_glm_pool(data_sets[[data_name]], f_cov, exposure)
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


# dist_plot <- bind_rows(
#     true_betas %>% mutate(method = "True"),
#     cc_betas %>% mutate(method = "Complete Case"),
#     imp_betas %>% mutate(method = "Multiple Imputation")
# ) |>
#     ggplot(aes(x = est, color = method)) +
#     geom_density(linewidth = 1) +
#     geom_vline(data = means, aes(xintercept = est, color = method), linewidth = 1, linetype = 2) +
#     labs(
#         title = paste0("Density of ", exposure, " coefficient for ", outcome),
#         subtitle = paste0(
#             format(iterations, big.mark = ","), " samples of size ", format(sample_size, big.mark = ","), " with average missing ", exposure, " %: ",
#             miss_dat |>
#                 map_dbl(\(x) mean(is.na(x[[exposure]]))) |>
#                 mean() |>
#                 (\(x) trimws(format(round(x * 100, 1), nsmall = 1)))()
#         ),
#         x = paste0(exposure, " coefficient"),
#         y = "Density"
#     ) +
#     ms::scale_color_ms() +
#     ms::theme_ms() +
#     theme(
#         legend.title = element_blank()
#     )

# ggsave(paste0("results/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, ".pdf"), plot = dist_plot, width = 8, height = 6, device = cairo_pdf)

# library(gridExtra)
# tt <- ttheme_default(
#     colhead = list(fg_params = list(parse = TRUE)),
#     base_size = 10,
#     padding = unit(c(2, 4), "mm")
# )
# tbl <- tableGrob(diag_tab, rows = NULL, theme = tt)

# cairo_pdf(paste0("results/", outcome, "_", exposure, "_n", sample_size, "_i", iterations, "_w_table.pdf"), width = 8, height = 6)
# grid.arrange(dist_plot, tbl,
#     nrow = 2, heights = c(2, 0.5)
# )
# dev.off()