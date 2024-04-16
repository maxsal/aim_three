# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr, glue
)

outcome_phe <- "EM_202.2"

plan(multicore, workers = 10)

# data -------------------------------------------------------------------------
analytic_sub <- qread(glue("data/private/analytic_sub_{outcome_phe}.qs"))
pim          <- qread(glue("data/private/mgi/20220822/MGI_PIM0X_20220822_{outcome_phe}.qs"))
ex_prs       <- qread("data/private/ex_prs_processed.qs")

# keep these pim vars
pim_sub_vars <- c(
    "id" = "id",
    "case" = "case"
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
        case # outcome
    )
d_id <- d |> select(id)
d <- d |> select(-id)

left_join(
    cbind(d_id, d) |> select(id, age, female, bmi, smoke, case),
    prs |> select(id, bmi_prs),
    by = "id"
) |>
    select(-id) |>
    cor()

# true beta
truth <- glm(bmi ~ case, data = d)$coefficients[["case"]]

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
