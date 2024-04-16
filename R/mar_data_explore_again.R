# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr, glue
)

outcome_phe <- "EM_202.2"

plan(multicore, workers = 10)

# data -------------------------------------------------------------------------
# analytic_sub <- qread(glue("data/private/analytic_sub_{outcome_phe}.qs"))
analytic_sub <- qread(glue("data/private/analytic_sub.qs"))
# pim          <- qread(glue("data/private/mgi/20220822/MGI_PIM0X_20220822_{outcome_phe}.qs"))
pim <- qread(glue("data/private/mgi/20220822/MGI_PIM0X_20220822.qs"))
ex_prs       <- qread("data/private/ex_prs_processed.qs")

# keep these pim vars
pim_sub_vars <- c(
    "id" = "id",
    "t2d" = outcome_phe
)

cond <- pim |>
    select(any_of(pim_sub_vars))

# keep these prs vars
prs_sub_vars <- c(
    "id" = "id",
    "bmi_prs" = "bmi",
    "t2d_prs" = "t2d"
)

prs <- ex_prs |>
    mutate(
        bmi_prs = scale(bmi)[, 1],
        t2d_prs = scale(t2d)[, 1]
    )

# merge pim vars into analytic_sub
merged <- analytic_sub |> left_join(cond, by = "id")

# keep limited data
d <- merged |>
    mutate(nhw = as.numeric(race_eth == "NHW")) |>
    select(
        id,
        age, female,
        bmi = bmi_med, smoke = smoke_ever,
        glucose = glucose_med,
        t2d
    ) |>
    mutate(
        obese = as.numeric(bmi >= 30)
    )
d_id <- d |> select(id)
d <- d |> select(-id)

cor_plot <- ggcorrplot(
    cor(d),
    lab = TRUE,
    show.diag = FALSE,
    colors = c("#0072B2", "white", "#D55E00")
)

ggsave(
    plot = cor_plot,
    filename = "bin/cor_plot.pdf",
    device = cairo_pdf,
    width = 5, height = 5
)

covs <- c("age", "female", "smoke")

# GLUCOSE ~ BMI
outcome  <- "glucose"
exposure <- "bmi"
var_miss <- "glucose"
family <- "gaussian"
f <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(c(exposure, covs), collapse = " + "))

d_copy <- d |> select(all_of(c(outcome, exposure, covs)))

un_truth <- glm(f, data = d_copy, family = family)$coefficients[[exposure]]
adj_truth <- glm(f_cov, data = d_copy, family = family)$coefficients[[exposure]]

true_betas <- tibble(
    betas = c(un_truth, adj_truth),
    cov = c("Unadjusted", "Covariate-adjusted")
)

make_mar_data <- function(data, var, covs, wgts = NULL, int_wgt = -1.5, print = TRUE) {
    data_matrix <- as.matrix(cbind(1, scale(data |> dplyr::select(any_of(covs)))))
    if (is.null(wgts)) {
        wgts <- rep(1, length(covs))
    } else if (length(wgts) != length(covs)) {
        stop("Length of weights (`wgts`) must match length of covariates (`covs`)")
    }
    score <- data_matrix %*% c(int_wgt, wgts)
    miss_prob <- plogis(score)
    miss_ind <- rbinom(nrow(data), 1, miss_prob)
    if (print) message(paste0("Miss%: ", round(mean(miss_ind) * 100, 1)))
    data[[var]][miss_ind == 1] <- NA
    return(list(
        data = data,
        prob = miss_prob
    ))
}

mar_data <- future_map(
    1:100,
    \(x) {
        make_mar_data(
            data = d_copy,
            var = var_miss,
            covs = c(covs, setdiff(c(exposure, outcome), var_miss)),
            wgts = c(3, 0.5, -1, 1),
            print = FALSE
        )
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

miss_prob_plot <- tibble(
    "prob" = mar_data[[1]][["prob"]]
) |>
    ggplot(aes(x = prob)) +
    geom_histogram(aes(y = ..ncount..)) +
    geom_density(aes(y = ..scaled..), color = "#0072B2", linewidth = 1) +
    labs(
        title = "Missingness probability distribution",
        x = "Probability",
        y = "Density"
    ) +
    theme_ms()
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_miss_prob_dist.pdf"),
    plot = miss_prob_plot,
    device = cairo_pdf,
    width = 5, height = 5
)

mar_cc_un_betas <- future_map_dbl(
    mar_data,
    \(x) {
        glm(
            f,
            data = x[["data"]],
            family = family
        )$coefficients[[exposure]]
    },
    .progress = TRUE
)

mar_cc_adj_betas <- future_map_dbl(
    mar_data,
    \(x) {
        glm(
            f_cov,
            data = x[["data"]],
            family = family
        )$coefficients[[exposure]]
    },
    .progress = TRUE
)

mar_imp_data <- future_map(
    mar_data,
    \(x) {
        mice(x[["data"]], m = 10, printFlag = FALSE)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

make_imp_dist_plot <- function(imp_data, var) {
    if (!is.numeric(imp_data[["data"]][[var]])) stop("Variable must be numeric")
    miss_index <- which(is.na(imp_data[["data"]][[var]]))
    no_miss <- tibble(
        no_miss_vals = imp_data[["data"]] |>
            slice((1:nrow(imp_data[["data"]]))[-miss_index]) |>
            pull(var)
    )
    out_plot <- no_miss |>
        ggplot(aes(x = no_miss_vals)) +
        geom_line(stat = "density", color = "#0072B2", linewidth = 1) +
        labs(
            title = paste0("Imputed ", var, " distribution"),
            x = var,
            y = "Density"
        ) +
        theme_ms()
    for (i in 1:imp_data[["m"]]) {
        imp_miss <- tibble(
            imp_miss_vals = complete(imp_data, i) |>
            slice(miss_index) |>
            pull(var)
        )
        out_plot <- out_plot +
            geom_line(
                data = imp_miss,
                aes(x = imp_miss_vals),
                stat = "density",
                color = "#D55E00", linewidth = 1, alpha = 0.5)
    }
    out_plot
}


these_iter <- sample(1:100, 6, replace = FALSE)
imp_dist_plots <- map(
    these_iter,
    \(i) {
        make_imp_dist_plot(mar_imp_data[[i]], var_miss)
    }
)
imp_patch_plot <- wrap_plots(imp_dist_plots, ncol = 3, nrow = 2)
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_imp_dist.pdf"),
    plot = imp_patch_plot,
    device = cairo_pdf,
    width = 10, height = 8
)

make_imp_scat_plot <- function(imp_data, miss_var, y_var, obs = 250) {
    if (!is.numeric(imp_data[["data"]][[miss_var]])) stop("Variable must be numeric")
    miss_index <- which(is.na(imp_data[["data"]][[miss_var]]))
    no_miss <- imp_data[["data"]] |>
        slice((1:nrow(imp_data[["data"]]))[-miss_index]) |>
        slice_sample(n = obs)
    out_plot <- no_miss |>
        ggplot(aes(x = .data[[miss_var]], y = .data[[y_var]])) +
        geom_point(color = "#0072B2", alpha = 0.25) +
        labs(
            title = paste0("Imputed vs. observed"),
            x = miss_var,
            y = y_var
        ) +
        theme_ms()
    out_plots <- list() 
    for (i in 1:imp_data[["m"]]) {
        imp_miss <- complete(imp_data, i) |>
            slice(miss_index) |>
            slice_sample(n = obs)
        out_plots[[i]] <- out_plot +
            geom_point(
                data = imp_miss,
                aes(x = .data[[miss_var]], y = .data[[y_var]]),
                color = "#D55E00", alpha = 0.25
            ) +
            geom_smooth(
                data = no_miss,
                aes(x = .data[[miss_var]], y = .data[[y_var]]),
                method = "loess", se = FALSE, formula = y ~ x,
                color = "#0072B2", linewidth = 1, span = 0.75
            ) +
            geom_smooth(
                data = imp_miss,
                aes(x = .data[[miss_var]], y = .data[[y_var]]),
                method = "loess", se = FALSE, formula = y ~ x,
                color = "#D55E00", linewidth = 1, span = 0.75
            )
    }
    patchwork::wrap_plots(plotlist = out_plots, ncol = 5, nrow = 2)
}
imp_scat_plot <- make_imp_scat_plot(mar_imp_data[[3]], var_miss, outcome)
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_imp_scat.png"),
    plot = imp_scat_plot,
    width = 12, height = 8, units = "in",
    dpi = 240
)

mar_imp_un_betas <- future_map_dbl(
    mar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    f,
                    data = complete(x, i),
                    family = family
                )$coefficients[[exposure]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mar_imp_adj_betas <- future_map_dbl(
    mar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    f_cov,
                    data = complete(x, i),
                    family = family
                )$coefficients[[exposure]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

beta_dist_plot <- bind_rows(
    tibble(betas = mar_cc_un_betas, type = "Complete case", cov = "Unadjusted"),
    tibble(betas = mar_cc_adj_betas, type = "Complete case", cov = "Covariate-adjusted"),
    tibble(betas = mar_imp_un_betas, type = "Imputed", cov = "Unadjusted"),
    tibble(betas = mar_imp_adj_betas, type = "Imputed", cov = "Covariate-adjusted")
) |>
    ggplot(aes(x = betas, color = type, linetype = type)) +
    geom_vline(data = true_betas, aes(xintercept = betas), color = "black", linewidth = 1) +
    geom_density(linewidth = 1) +
    labs(
        title = paste0("BMI (cont) beta coefficient for glucose (cont) (", f, ")"),
        subtitle = paste0(
            "100 datasets; average ", var_miss, " missing: ",
            round(mean(map_dbl(mar_data, \(x) mean(is.na(x[["data"]][[exposure]]))) * 100, 1)),
            "%"
            ),
        x = "Beta coefficient",
        y = "Density"
    ) +
    scale_color_ms() +
    facet_wrap(~factor(cov, levels = c("Unadjusted", "Covariate-adjusted")), ncol = 1) +
    theme_ms()
ggsave(
    paste0("bin/",outcome, "_", exposure, "_", var_miss, "m_beta_dist.pdf"),
    plot = beta_dist_plot,
    device = cairo_pdf,
    width = 8, height = 6
)

outcome <- "glucose"
exposure <- "bmi"
var_miss <- "glucose"
family <- "gaussian"
f <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(c(exposure, covs), collapse = " + "))

d_copy <- d |> select(all_of(c(outcome, exposure, covs)))

un_truth <- glm(f, data = d_copy, family = family)$coefficients[[exposure]]
adj_truth <- glm(f_cov, data = d_copy, family = family)$coefficients[[exposure]]

true_betas <- tibble(
    betas = c(un_truth, adj_truth),
    cov = c("Unadjusted", "Covariate-adjusted")
)

mar_data <- future_map(
    1:100,
    \(x) {
        make_mar_data(
            data = d_copy,
            var = var_miss,
            covs = c(covs, setdiff(c(exposure, outcome), var_miss)),
            wgts = c(3, 0.5, -1, 1),
            print = FALSE
        )
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

miss_prob_plot <- tibble(
    "prob" = mar_data[[1]][["prob"]]
) |>
    ggplot(aes(x = prob)) +
    geom_histogram(aes(y = ..ncount..)) +
    geom_density(aes(y = ..scaled..), color = "#0072B2", linewidth = 1) +
    labs(
        title = "Missingness probability distribution",
        x = "Probability",
        y = "Density"
    ) +
    theme_ms()
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_miss_prob_dist.pdf"),
    plot = miss_prob_plot,
    device = cairo_pdf,
    width = 5, height = 5
)

mar_cc_un_betas <- future_map_dbl(
    mar_data,
    \(x) {
        glm(
            f,
            data = x[["data"]],
            family = family
        )$coefficients[[exposure]]
    },
    .progress = TRUE
)

mar_cc_adj_betas <- future_map_dbl(
    mar_data,
    \(x) {
        glm(
            f_cov,
            data = x[["data"]],
            family = family
        )$coefficients[[exposure]]
    },
    .progress = TRUE
)

mar_imp_data <- future_map(
    mar_data,
    \(x) {
        mice(x[["data"]], m = 10, printFlag = FALSE)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

these_iter <- sample(1:100, 6, replace = FALSE)
imp_dist_plots <- map(
    these_iter,
    \(i) {
        make_imp_dist_plot(mar_imp_data[[i]], var_miss)
    }
)
imp_patch_plot <- wrap_plots(imp_dist_plots, ncol = 3, nrow = 2)
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_imp_dist.pdf"),
    plot = imp_patch_plot,
    device = cairo_pdf,
    width = 10, height = 8
)

imp_scat_plot <- make_imp_scat_plot(mar_imp_data[[3]], var_miss, setdiff(c(exposure, outcome), var_miss))
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_imp_scat.png"),
    plot = imp_scat_plot,
    width = 12, height = 8, units = "in",
    dpi = 240
)

mar_imp_un_betas <- future_map_dbl(
    mar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    f,
                    data = complete(x, i),
                    family = family
                )$coefficients[[exposure]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mar_imp_adj_betas <- future_map_dbl(
    mar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    f_cov,
                    data = complete(x, i),
                    family = family
                )$coefficients[[exposure]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

beta_dist_plot <- bind_rows(
    tibble(betas = mar_cc_un_betas, type = "Complete case", cov = "Unadjusted"),
    tibble(betas = mar_cc_adj_betas, type = "Complete case", cov = "Covariate-adjusted"),
    tibble(betas = mar_imp_un_betas, type = "Imputed", cov = "Unadjusted"),
    tibble(betas = mar_imp_adj_betas, type = "Imputed", cov = "Covariate-adjusted")
) |>
    ggplot(aes(x = betas, color = type, linetype = type)) +
    geom_vline(data = true_betas, aes(xintercept = betas), color = "black", linewidth = 1) +
    geom_density(linewidth = 1) +
    labs(
        title = paste0("BMI (cont) beta coefficient for glucose (cont) (", f, ")"),
        subtitle = paste0(
            "100 datasets; average ", var_miss, " missing: ",
            round(mean(map_dbl(mar_data, \(x) mean(is.na(x[["data"]][[exposure]]))) * 100, 1)),
            "%"
        ),
        x = "Beta coefficient",
        y = "Density"
    ) +
    scale_color_ms() +
    facet_wrap(~ factor(cov, levels = c("Unadjusted", "Covariate-adjusted")), ncol = 1) +
    theme_ms()
ggsave(
    paste0("bin/", outcome, "_", exposure, "_", var_miss, "m_beta_dist.pdf"),
    plot = beta_dist_plot,
    device = cairo_pdf,
    width = 8, height = 6
)

# introduce missingness --------------------------------------------------------
# MCAR
make_mcar_data <- function(data, var, prop = 0.3, outcome = "case") {
    miss_prop <- rbinom(nrow(data), 1, prop)
    data[[var]][miss_prop == 1] <- NA
    data
}

mcar_data <- future_map(
    1:100,
    \(x) {
        make_mcar_data(data = d, var = "bmi", prop = 0.3)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mcar_cc_betas <- future_map_dbl(
    mcar_data,
    \(x) {
        glm(
            bmi ~ case,
            data = x
        )$ coefficients[["case"]]
        # logistf::logistf(
        #     case ~ bmi,
        #     data = x,
        #     pl = FALSE,
        #     control = logistf::logistf.control(
        #         maxit = 1000,
        #         maxstep = 0.5
        #     )
        # )$coefficients[["bmi"]]
    },
    .progress = TRUE
)

mcar_imp_data <- future_map(
    mcar_data,
    \(x) {
        mice(x, m = 10, printFlag = FALSE)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mcar_imp_betas <- future_map_dbl(
    mcar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    bmi ~ case,
                    data = complete(x, i)
                )$coefficients[["case"]]
                # logistf::logistf(
                #     case ~ bmi,
                #     data = complete(x, i),
                #     pl = FALSE,
                #     control = logistf::logistf.control(
                #         maxit = 1000,
                #         maxstep = 0.5
                #     )
                # )$coefficients[["bmi"]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mcar_imp_bmi_prs_data <- future_map(
    mcar_data,
    \(x) {
        cbind(d_id, x) |>
            left_join(prs |> dplyr::select(id, bmi_prs), by = "id") |>
            select(-id) |>
            mice(m = 10, printFlag = FALSE)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mcar_imp_bmi_prs_betas <- future_map_dbl(
    mcar_imp_bmi_prs_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    bmi ~ case,
                    data = complete(x, i)
                )$coefficients[["case"]]
                # logistf::logistf(
                #     case ~ bmi,
                #     data = complete(x, i),
                #     pl = FALSE,
                #     control = logistf::logistf.control(
                #         maxit = 1000,
                #         maxstep = 0.5
                #     )
                # )$coefficients[["bmi"]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

# MAR
# step 1. create the datasets
# step 2. create the scores
# step 3. create the missing probabilities
# step 4. create the missing indicators
# step 5. estimate complete case beta coefficients
# step 6. impute the missing values
# step 7. estimate the imputed beta coefficients
make_mar_data <- function(data, var, covs, wgts = NULL, int_wgt = -1.5, print = TRUE) {
    data_matrix <- as.matrix(cbind(1, scale(data |> dplyr::select(any_of(covs)))))
    if (is.null(wgts)) {
        wgts <- rep(1, length(covs))
    } else if (length(wgts) != length(covs)) {
        stop("Length of weights (`wgts`) must match length of covariates (`covs`)")
    }
    score <- data_matrix %*% c(int_wgt, wgts)
    miss_prob <- plogis(score)
    miss_ind <- rbinom(nrow(data), 1, miss_prob)
    if (print) message(paste0("Miss%: ", round(mean(miss_ind) * 100, 1)))
    data[[var]][miss_ind == 1] <- NA
    return(data)
}

mar_data <- future_map(
    1:100,
    \(x) {
        make_mar_data(
            data = d,
            var = "bmi",
            covs = c("age", "female", "smoke", "case"),
            wgts = c(3, 0.5, -1, 1),
            print = FALSE
        )
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mar_cc_betas <- future_map_dbl(
    mar_data,
    \(x) {
        glm(
            bmi ~ case,
            data = x
        )$coefficients[["case"]]
        # logistf::logistf(
        #     case ~ bmi,
        #     data = x,
        #     pl = FALSE,
        #     control = logistf::logistf.control(
        #         maxit = 1000,
        #         maxstep = 0.5)
        #     )$coefficients[["bmi"]]
    },
    .progress = TRUE
)

mar_imp_data <- future_map(
    mar_data,
    \(x) {
        mice(x, m = 10, printFlag = FALSE, method = "lasso.norm")
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mar_imp_betas <- future_map_dbl(
    mar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    bmi ~ case,
                    data = complete(x, i)
                )$coefficients[["case"]]
                # logistf::logistf(
                #     case ~ bmi,
                #     data = complete(x, i),
                #     pl = FALSE,
                #     control = logistf::logistf.control(
                #         maxit = 1000,
                #         maxstep = 0.5
                #     )
                # )$coefficients[["bmi"]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mar_imp_bmi_prs_data <- future_map(
    mar_data,
    \(x) {
        cbind(d_id, x) |>
            left_join(prs |> dplyr::select(id, bmi_prs), by = "id") |>
            select(-id) |>
            mice(m = 10, printFlag = FALSE, method = "lasso.norm")
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mar_imp_bmi_prs_betas <- future_map_dbl(
    mar_imp_bmi_prs_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    bmi ~ case,
                    data = complete(x, i)
                )$coefficients[["case"]]
                # logistf::logistf(
                #     case ~ bmi,
                #     data = complete(x, i),
                #     pl = FALSE,
                #     control = logistf::logistf.control(
                #         maxit = 1000,
                #         maxstep = 0.5
                #     )
                # )$coefficients[["bmi"]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

# mar_imp_bmi_t2d_prs_data <- future_map(
#     mar_data,
#     \(x) {
#         cbind(d_id, x) |>
#             left_join(prs |> dplyr::select(id, bmi_prs, t2d_prs), by = "id") |>
#             select(-id) |>
#             mice(m = 10, printFlag = FALSE)
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )

# mar_imp_bmi_t2d_prs_betas <- future_map_dbl(
#     mar_imp_bmi_t2d_prs_data,
#     \(x) {
#         betas <- map_dbl(
#             1:10,
#             \(i) {
#                 logistf::logistf(
#                     case ~ bmi,
#                     data = complete(x, i),
#                     pl = FALSE,
#                     control = logistf::logistf.control(
#                         maxit = 1000,
#                         maxstep = 0.5
#                     )
#                 )$coefficients[["bmi"]]
#             }
#         )
#         mean(betas)
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )

# MNAR
make_mnar_data <- function(data, var, covs, wgts = NULL, int_wgt = -1.5, var_wgt = 2, print = TRUE) {
    data_matrix <- as.matrix(cbind(1, scale(data |> dplyr::select(any_of(c(covs, var))))))
    if (is.null(wgts)) {
        wgts <- rep(1, length(covs))
    } else if (length(wgts) != length(covs)) {
        stop("Length of weights (`wgts`) must match length of covariates (`covs`)")
    }
    score <- data_matrix %*% c(int_wgt, wgts, var_wgt)
    miss_prob <- plogis(score)
    miss_ind <- rbinom(nrow(data), 1, miss_prob)
    if (print) message(paste0("Miss%: ", round(mean(miss_ind) * 100, 1)))
    data[[var]][miss_ind == 1] <- NA
    return(data)
}

mnar_data <- future_map(
    1:100,
    \(x) {
        make_mnar_data(
            data = d,
            var = "bmi",
            covs = c("age", "female", "smoke", "case"),
            wgts = c(3, 0.5, -1, 1),
            print = FALSE
        )
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mnar_cc_betas <- future_map_dbl(
    mnar_data,
    \(x) {
        glm(
            bmi ~ case,
            data = x
        )$coefficients[["case"]]
        # logistf::logistf(
        #     case ~ bmi,
        #     data = x,
        #     pl = FALSE,
        #     control = logistf::logistf.control(
        #         maxit = 1000,
        #         maxstep = 0.5
        #     )
        # )$coefficients[["bmi"]]
    },
    .progress = TRUE
)

mnar_imp_data <- future_map(
    mnar_data,
    \(x) {
        mice(x, m = 10, printFlag = FALSE)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mnar_imp_betas <- future_map_dbl(
    mnar_imp_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    bmi ~ case,
                    data = complete(x, i)
                )$coefficients[["case"]]
                # logistf::logistf(
                #     case ~ bmi,
                #     data = complete(x, i),
                #     pl = FALSE,
                #     control = logistf::logistf.control(
                #         maxit = 1000,
                #         maxstep = 0.5
                #     )
                # )$coefficients[["bmi"]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mnar_imp_bmi_prs_data <- future_map(
    mnar_data,
    \(x) {
        cbind(d_id, x) |>
            left_join(prs |> dplyr::select(id, bmi_prs), by = "id") |>
            select(-id) |>
            mice(m = 10, printFlag = FALSE)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

mnar_imp_bmi_prs_betas <- future_map_dbl(
    mnar_imp_bmi_prs_data,
    \(x) {
        betas <- map_dbl(
            1:10,
            \(i) {
                glm(
                    bmi ~ case,
                    data = complete(x, i)
                )$coefficients[["case"]]
                # logistf::logistf(
                #     case ~ bmi,
                #     data = complete(x, i),
                #     pl = FALSE,
                #     control = logistf::logistf.control(
                #         maxit = 1000,
                #         maxstep = 0.5
                #     )
                # )$coefficients[["bmi"]]
            }
        )
        mean(betas)
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

# mnar_imp_bmi_t2d_prs_data <- future_map(
#     mnar_data,
#     \(x) {
#         cbind(d_id, x) |>
#             left_join(prs |> dplyr::select(id, bmi_prs, t2d_prs), by = "id") |>
#             select(-id) |>
#             mice(m = 10, printFlag = FALSE)
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )

# mnar_imp_bmi_t2d_prs_betas <- future_map_dbl(
#     mnar_imp_bmi_t2d_prs_data,
#     \(x) {
#         betas <- map_dbl(
#             1:10,
#             \(i) {
#                 logistf::logistf(
#                     case ~ bmi,
#                     data = complete(x, i),
#                     pl = FALSE,
#                     control = logistf::logistf.control(
#                         maxit = 1000,
#                         maxstep = 0.5
#                     )
#                 )$coefficients[["bmi"]]
#             }
#         )
#         mean(betas)
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )

# quick summary plot -----------------------------------------------------------
bind_rows(
    tibble(
        betas = mcar_cc_betas,
        type = "Complete case",
        mech = "MCAR"
    ),
    tibble(
        betas = mcar_imp_betas,
        type = "OD-imputed",
        mech = "MCAR"
    ),
    tibble(
        betas = mcar_imp_bmi_prs_betas,
        type = "BMIprs-imputed",
        mech = "MCAR"
    ),
    tibble(
        betas = mar_cc_betas,
        type = "Complete case",
        mech = "MAR"
    ),
    tibble(
        betas = mar_imp_betas,
        type = "OD-imputed",
        mech = "MAR"
    ),
    tibble(
        betas = mar_imp_bmi_prs_betas,
        type = "BMIprs-imputed",
        mech = "MAR"
    ),
    # tibble(
    #     betas = mar_imp_bmi_t2d_prs_betas,
    #     type = "BMIT2Dprs-imputed",
    #     mech = "MAR"
    # ),
    tibble(
        betas = mnar_cc_betas,
        type = "Complete case",
        mech = "MNAR"
    ),
    tibble(
        betas = mnar_imp_betas,
        type = "OD-imputed",
        mech = "MNAR"
    ),
    tibble(
        betas = mnar_imp_bmi_prs_betas,
        type = "BMIprs-imputed",
        mech = "MNAR"
    )
    # tibble(
    #     betas = mnar_imp_bmi_t2d_prs_betas,
    #     type = "BMIT2Dprs-imputed",
    #     mech = "MNAR"
    # )
) |>
    mutate(
        `Missingness mechanism` = factor(mech, levels = c("MCAR", "MAR", "MNAR")),
        `Analysis` = factor(type, levels = c("Complete case", "OD-imputed", "BMIprs-imputed", "BMIT2Dprs-imputed"))
    ) |>
    ggplot(aes(x = betas, color = `Missingness mechanism`, linetype = `Analysis`)) +
    geom_vline(xintercept = truth, linewidth = 1) +
    geom_density(linewidth = 1) +
    ms::scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.box = "vertical"
    )

ggsave(
    "bin/missing_dist.pdf",
    device = cairo_pdf,
    width = 10, height = 4
)

bind_rows(
    tibble(
        betas = mar_cc_betas,
        type = "Complete case",
        mech = "MAR"
    ),
    tibble(
        betas = mar_imp_betas,
        type = "OD-imputed",
        mech = "MAR"
    ),
    tibble(
        betas = mar_imp_bmi_prs_betas,
        type = "BMIprs-imputed",
        mech = "MAR"
    ),
    tibble(
        betas = mar_imp_bmi_t2d_prs_betas,
        type = "BMIT2Dprs-imputed",
        mech = "MAR"
    )
) |>
    mutate(
        `Missingness mechanism` = factor(mech, levels = c("MCAR", "MAR", "MNAR")),
        `Analysis` = factor(type, levels = c("Complete case", "OD-imputed", "BMIprs-imputed", "BMIT2Dprs-imputed"))
    ) |>
    ggplot(aes(x = betas, color = `Missingness mechanism`, linetype = `Analysis`)) +
    geom_vline(xintercept = truth, linewidth = 1) +
    geom_density(linewidth = 1) +
    scale_color_manual(
        values = c(
            "MCAR" = "#E69F00",
            "MAR" = "#56B4E9",
            "MNAR" = "#009E73"
        )
    ) +
    ms::theme_ms() +
    theme(
        legend.box = "vertical"
    )

ggsave(
    "bin/mar_dist.pdf",
    device = cairo_pdf,
    width = 10, height = 4
)
