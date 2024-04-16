# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr
)

plan(multicore, workers = 10)

# data -------------------------------------------------------------------------
analytic_sub <- qread("data/private/analytic_sub.qs")
pim <- qread("data/private/mgi/20220822/MGI_PIM0X_20220822.qs")

# keep these pim vars
pim_sub_vars <- c(
    "id" = "id",
    "t2d" = "EM_202.2",
    "hyper" = "CV_401",
    "cor_ath" = "CV_404.2",
    "crc" = "CA_101.41",
    "liv_can" = "CA_101.6",
    "pan_can" = "CA_101.8",
    "kid_can" = "CA_108.4"

)

cond <- pim |>
    select(any_of(pim_sub_vars))

# merge pim vars into analytic_sub
merged <- analytic_sub |> left_join(cond, by = "id")

# keep limited data
limited <- merged |>
    select(
        age, female, race_eth, cci = charlson_quan_weighted, # covariates
        height = height_cm_med, bmi = bmi_med, smoke = smoke_ever, glucose = glucose_med, alcohol = alcohol_ever, # exposures
        t2d, hyper, cor_ath # outcomes
    ) |>
    mutate(
        obese = ifelse(bmi >= 30, 1, 0)
    ) |> select(age, female, obese, cci, t2d)

expanded <- merged |>
    mutate(nhw = as.numeric(race_eth == "NHW")) |>
    select(
        age, female, nhw, cci = charlson_quan_weighted, # covariates
        height = height_cm_med, bmi = bmi_med, smoke = smoke_ever, glucose = glucose_med, alcohol = alcohol_ever, # exposures
        t2d, hyper, cor_ath # outcomes
    )

# explore ampute base patterns
base_pattern <- ampute(expanded)$patterns
base_pattern$age <- base_pattern$female <- base_pattern$nhw <-  base_pattern$cci <- base_pattern$t2d <- base_pattern$hyper <- base_pattern$cor_ath <- 1
keep <- apply(base_pattern, 1, \(x) length(unique(x[!is.na(x)])) != 1)
base_pattern[keep, ]

new_patterns <- rbind(
    base_pattern[keep, ],
    c(1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1),
    c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
    c(1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1)
)

result <- ampute(expanded, patterns = new_patterns)

patt <- rep(1, ncol(limited))
patt[which(names(limited) == "obese")] <- 0
wts <- rep(1, ncol(limited))
wts[which(names(limited) %in% c("obese", "t2d"))] <- 0

# map_dbl(mar_0.5$amp, \(x) sum(is.na(x))) |> enframe() |> arrange(desc(value))

# setup for ampute/impute function
outcome   <- "t2d"
predictor <- .predictor <- "bmi"
covs      <- c("age", "female", "nhw", "cci")
exposures <- c("height", "bmi", "smoke", "glucose", "alcohol")
exposures <- exposures[!exposures %in% predictor]
n_seeds   <- 200
n_imps    <- 10
miss_prop <- 0.5

# complete case
# unadjusted
un_adj_cc <- glm(paste0(outcome, "~", predictor), data = expanded, family = "binomial") |> aimTwo::betas_from_mod()
# covariate adjusted
cov_adj_cc <- glm(paste0(outcome, "~", paste0(c(predictor, covs), collapse = "+")), data = expanded, family = "binomial") |> aimTwo::betas_from_mod()
# full adjusted
full_adj_cc <- glm(paste0(outcome, "~", paste0(c(predictor, covs, exposures), collapse = "+")), data = expanded, family = "binomial") |> aimTwo::betas_from_mod()

vlines <- tibble(
    adjust = c("un_adj", "cov_adj", "full_adj"),
    beta = c(
        un_adj_cc |> filter(predictor == .predictor) |> pull(beta),
        cov_adj_cc |> filter(predictor == .predictor) |> pull(beta),
        full_adj_cc |> filter(predictor == .predictor) |> pull(beta)
    )
)

cascade_amp_imp <- function(
    data,
    outcome,
    predictor,
    covs,
    exposures,
    patterns = NULL,
    freq  = NULL,
    weights = NULL,
    miss_prop = 0.5,
    mech = "MAR",
    n_seeds = 10,
    n_imp = 10,
    .progress = TRUE
) {
    future_map_dfr(
        1:n_seeds,
        \(x) {
            set.seed(x)
            tmp <- ampute(data, patterns = patterns, freq = freq, weights = weights, prop = miss_prop, mech = mech)
            .predictor <- predictor

            # formula
            un_f <- paste0(outcome, "~", predictor)
            cov_f <- paste0(outcome, "~", paste0(c(predictor, covs), collapse = "+"))
            full_f <- paste0(outcome, "~", paste0(c(predictor, covs, exposures), collapse = "+"))


            ## complete case
            # unadjusted
            un_adj <- glm(un_f, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                filter(predictor == .predictor) |>
                mutate(seed = x, adjust = "un_adj", type = "ampute_cc")
            # covariate adjusted
            cov_adj <- glm(cov_f, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                filter(predictor == .predictor) |>
                mutate(seed = x, adjust = "cov_adj", type = "ampute_cc")
            # full adjusted
            full_adj <- glm(full_f, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                filter(predictor == .predictor) |>
                mutate(seed = x, adjust = "full_adj", type = "ampute_cc")

            ## imputed: observed cases
            # impute
            tmp_imp <- mice(tmp$amp, m = n_imp, print = FALSE)
            # unadjusted
            pool1 <- pool(with(tmp_imp, glm(formula(un_f), family = "binomial")))
            un_adj_imp <- as.data.table(pool1$pooled)[term == predictor, .(predictor = term, beta = estimate)] |>
                mutate(seed = x, adjust = "un_adj", type = "imputed_oc")
            # covariate adjusted
            pool2 <- pool(with(tmp_imp, glm(formula(cov_f), family = "binomial")))
            cov_adj_imp <- as.data.table(pool2$pooled)[term == predictor, .(predictor = term, beta = estimate)] |>
                mutate(seed = x, adjust = "cov_adj", type = "imputed_oc")
            # covariate adjusted
            pool3 <- pool(with(tmp_imp, glm(formula(full_f), family = "binomial")))
            full_adj_imp <- as.data.table(pool3$pooled)[term == predictor, .(predictor = term, beta = estimate)] |>
                mutate(seed = x, adjust = "full_adj", type = "imputed_oc")

            ## imputed: observed cases + predictor ExPRS
            ## imputed: observed cases + all ExPRS

            # combine
            bind_rows(
                un_adj, cov_adj, full_adj,
                un_adj_imp, cov_adj_imp, full_adj_imp
            )
        },
        .progress = .progress,
        .options = furrr_options(seed = TRUE)
    )
}


mcar <- cascade_amp_imp(
    data = expanded,
    outcome = "t2d",
    predictor = "bmi",
    covs = covs,
    exposures = exposures,
    patterns = new_patterns,
    miss_prop = miss_prop,
    mech = "MCAR",
    n_seeds = n_seeds,
    n_imp = n_imps,
    .progress = TRUE
)

mar <- cascade_amp_imp(
    data = expanded,
    outcome = "t2d",
    predictor = "bmi",
    covs = covs,
    exposures = exposures,
    patterns = new_patterns,
    miss_prop = miss_prop,
    mech = "MAR",
    n_seeds = n_seeds,
    n_imp = n_imps,
    .progress = TRUE
)

mnar <- cascade_amp_imp(
    data      = expanded,
    outcome   = "t2d",
    predictor = "bmi",
    covs      = covs,
    exposures = exposures,
    patterns = new_patterns,
    miss_prop = miss_prop,
    mech = "MNAR",
    n_seeds = n_seeds,
    n_imp = n_imps,
    .progress = TRUE
)

test_plot <- bind_rows(
    mcar |> mutate(mech = "MCAR"),
    mar |> mutate(mech = "MAR"),
    mnar |> mutate(mech = "MNAR")
) |>
    ggplot(aes(x = beta, color = mech, linetype = type)) +
    geom_vline(data = vlines, aes(xintercept = beta)) +
    # geom_histogram(bins = 30, alpha = 0.5, position = "dodge") +
    geom_density(linewidth = 1, alpha = 0.5) +
    facet_wrap(~adjust, ncol = 1) +
    scale_color_ms() +
    theme_ms()

ggsave(
    "bin/test_plot.pdf",
    plot = test_plot,
    width = 6, height = 6, device = cairo_pdf
)


#### Some old code pre ampute/impute function
mar <- map_dfr(
    1:200,
    \(x) {
        set.seed(x)
        tmp <- ampute(limited, prop = 0.5, mech = "MAR")
        bind_rows(
            glm(t2d ~ obese, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                mutate(seed = x, adjust = "unadjusted"),
            glm(t2d ~ obese + age + female + cci, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                filter(predictor == "obese") |>
                mutate(seed = x, adjust = "adjusted")
        )
    },
    .progress = TRUE
)

mnar <- map_dfr(
    1:200,
    \(x) {
        set.seed(x)
        tmp <- ampute(limited, prop = 0.5, mech = "MNAR")
        bind_rows(
            glm(t2d ~ obese, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                mutate(seed = x, adjust = "unadjusted"),
            glm(t2d ~ obese + age + female + cci, data = tmp$amp, family = "binomial") |>
                aimTwo::betas_from_mod() |>
                filter(predictor == "obese") |>
                mutate(seed = x, adjust = "adjusted")
        )
    },
    .progress = TRUE
)

vlines <- tibble(
    adjust = c("unadjusted", "adjusted"),
    beta = c(1.00235, 1.0531765)
)

tmp_plot <- bind_rows(
    mcar |> mutate(mech = "MCAR"),
    mar |> mutate(mech = "MAR"),
    mnar |> mutate(mech = "MNAR")
) |>
    ggplot(aes(x = beta, color = mech)) +
    geom_vline(data = vlines, aes(xintercept = beta)) +
    geom_density(linewidth = 1, alpha = 0.5) +
    facet_wrap(~adjust, ncol = 1) +
    scale_color_ms() +
    theme_ms()

ggsave(
    "bin/mcar_mar_mnar_plot.pdf",
    plot = tmp_plot,
    width = 6, height = 6, device = cairo_pdf
)

bind_rows(
    mcar |> mutate(mech = "MCAR"),
    mar |> mutate(mech = "MAR"),
    mnar |> mutate(mech = "MNAR")
) |>
    ggplot(aes(x = beta, fill = mech)) +
    geom_vline(xintercept = 1.00235) +
    geom_histogram(bins = 30, alpha = 0.5, position = "dodge") +
    # geom_density(linewidth = 1, alpha = 0.5) +
    scale_fill_ms() +
    theme_ms()

mnar_0.5 <- ampute(limited, patterns = patt, weights = wts, prop = 0.5, mech = "MNAR")
mnar_0.5$weights
# mnar_0.5$amp |>
#     map_dbl(\(x) sum(is.na(x))) |>
#     enframe() |>
#     arrange(desc(value))
glm(t2d ~ obese, data = mnar_0.5$amp, family = "binomial") |> summary()


map_dfr(
    c("t2d", "hyper", "cor_ath", "kid_can"),
    \(x) {
        glm(
            as.formula(paste0(x, " ~ obese")),
            data = mnar_0.5$amp,
            family = "binomial"
        ) |>
            aimTwo::betas_from_mod() |>
            mutate(
                outcome = x
            )
    }
) |> select(outcome, predictor, or, p_value)


r <- rbinom(length(limited$obese), size = 1, prob = c(.25, .5)[limited$obese + 1])
limited$obese[r == 1] <- NA



map_dfr(
    c("t2d", "hyper", "cor_ath", "kid_can"),
    \(x) {
        glm(
            as.formula(paste0(x, " ~ bmi_med")),
            data = limited,
            family = "binomial"
        ) |>
            aimTwo::betas_from_mod() |>
            mutate(
                outcome = x
            )
    }
) |> select(outcome, predictor, or, p_value)

map_dfr(
    c("t2d", "hyper", "cor_ath", "kid_can"),
    \(x) {
        glm(
            as.formula(paste0(x, " ~ bmi_med")),
            data = mnar_0.5$amp,
            family = "binomial"
        ) |>
            aimTwo::betas_from_mod() |>
            mutate(
                outcome = x
            )
    }
) |> select(outcome, predictor, or, p_value)

mcar_0.3 <- ampute(limited, patterns = patt, prop = 0.3, mech = "MCAR")
mcar_0.5 <- ampute(limited, patterns = patt, prop = 0.5, mech = "MCAR")

mcar <- function(data, miss_var, miss_prop) {
    tmp <- data
    num_missing <- round(nrow(tmp) * miss_prop)
    miss_idx <- sample(1:nrow(tmp), num_missing)
    tmp[[miss_var]][miss_idx] <- NA
    return(tmp)
}
limited
mcar(limited, "bmi_med", miss_prop = 0.5)

glm(crc ~ bmi_med, data = taco, family = "binomial") |> summary()
glm(crc ~ bmi_med, data = mcar_0.5$amp, family = "binomial") |> summary()
glm(crc ~ bmi_med, data = mar_0.5$amp, family = "binomial") |> summary()
glm(crc ~ bmi_med, data = mnar_0.5$amp, family = "binomial") |> summary()


ampute(analytic_sub, prop = 0.3, mech = "MCAR")
ampute(analytic_sub, prop = 0.5, mech = "MCAR")

tmp <- ampute(analytic_sub |>
    select(age, female, bmi_med, cancerx), prop = 0.5, mech = "MNAR")
tmp |> names()

tmp$patterns |> dim()

map(
    seq_along(mechs),
    \(m) {
        map(
            seq_along(props),
            \(p) {
                ampute(analytic_sub, prop = props[p], mech = mechs[m])
            }
        )
    }
)


################################################################################
### OLD - FOR REFERENCE ########################################################
################################################################################

# mcar <- future_map_dfr(
#     1:10,
#     \(x) {
#         set.seed(x)
#         tmp <- ampute(expanded, prop = 0.5, mech = "MCAR")

#         # formula
#         un_f   <- paste0(outcome, "~", predictor)
#         cov_f  <- paste0(outcome, "~", paste0(c(predictor, covs), collapse = "+"))
#         full_f <- paste0(outcome, "~", paste0(c(predictor, covs, exposures), collapse = "+"))


#         ## complete case
#         # unadjusted
#         un_adj <- glm(un_f, data = tmp$amp, family = "binomial") |>
#             aimTwo::betas_from_mod() |>
#             filter(predictor == .predictor) |>
#             mutate(seed = x, adjust = "un_adj", type = "ampute_cc")
#         # covariate adjusted
#         cov_adj <- glm(cov_f, data = tmp$amp, family = "binomial") |>
#             aimTwo::betas_from_mod() |>
#             filter(predictor == .predictor) |>
#             mutate(seed = x, adjust = "cov_adj", type = "ampute_cc")
#         # full adjusted
#         full_adj <- glm(full_f, data = tmp$amp, family = "binomial") |>
#             aimTwo::betas_from_mod() |>
#             filter(predictor == .predictor) |>
#             mutate(seed = x, adjust = "full_adj", type = "ampute_cc")

#         ## imputed
#         # impute
#         tmp_imp <- mice(tmp$amp, m = 10, print = FALSE)
#         # unadjusted
#         pool1 <- pool(with(tmp_imp, glm(formula(un_f), family = "binomial")))
#         un_adj_imp <- as.data.table(pool1$pooled)[term == predictor, .(predictor = term, beta = estimate)] |>
#             mutate(seed = x, adjust = "un_adj", type = "imputed")
#         # covariate adjusted
#         pool2 <- pool(with(tmp_imp, glm(formula(cov_f), family = "binomial")))
#         cov_adj_imp <- as.data.table(pool2$pooled)[term == predictor, .(predictor = term, beta = estimate)] |>
#             mutate(seed = x, adjust = "cov_adj", type = "imputed")
#         # covariate adjusted
#         pool3 <- pool(with(tmp_imp, glm(formula(full_f), family = "binomial")))
#         full_adj_imp <- as.data.table(pool3$pooled)[term == predictor, .(predictor = term, beta = estimate)] |>
#             mutate(seed = x, adjust = "cov_adj", type = "imputed")

#         # combine
#         bind_rows(
#             un_adj, cov_adj, full_adj,
#             un_adj_imp, cov_adj_imp, full_adj_imp
#         )
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )

# mcar

# mar <- future_map_dfr(
#     1:100,
#     \(x) {
#         set.seed(x)
#         tmp <- ampute(limited, prop = 0.5, mech = "MAR")
#         tmp_imp <- mice(tmp$amp, m = 10, print = FALSE)
#         pool1 <- pool(with(tmp_imp, glm(t2d ~ obese, family = "binomial")))
#         pool2 <- pool(with(tmp_imp, glm(t2d ~ obese + age + female + cci, family = "binomial")))
#         as.data.table(pool1$pooled)[term == "obese", ]

#         bind_rows(
#             glm(t2d ~ obese, data = tmp$amp, family = "binomial") |>
#                 aimTwo::betas_from_mod() |>
#                 mutate(seed = x, adjust = "unadjusted", type = "complete"),
#             glm(t2d ~ obese + age + female + cci, data = tmp$amp, family = "binomial") |>
#                 aimTwo::betas_from_mod() |>
#                 filter(predictor == "obese") |>
#                 mutate(seed = x, adjust = "adjusted", type = "complete"),
#             as.data.table(pool1$pooled)[term == "obese", .(predictor = term, beta = estimate)] |>
#                 mutate(seed = x, adjust = "unadjusted", type = "imputed"),
#             as.data.table(pool2$pooled)[term == "obese", .(predictor = term, beta = estimate)] |>
#                 mutate(seed = x, adjust = "adjusted", type = "imputed")
#         )
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )

# mnar <- future_map_dfr(
#     1:100,
#     \(x) {
#         set.seed(x)
#         tmp <- ampute(limited, prop = 0.5, mech = "MNAR")
#         tmp_imp <- mice(tmp$amp, m = 10, print = FALSE)
#         pool1 <- pool(with(tmp_imp, glm(t2d ~ obese, family = "binomial")))
#         pool2 <- pool(with(tmp_imp, glm(t2d ~ obese + age + female + cci, family = "binomial")))
#         as.data.table(pool1$pooled)[term == "obese", ]

#         bind_rows(
#             glm(t2d ~ obese, data = tmp$amp, family = "binomial") |>
#                 aimTwo::betas_from_mod() |>
#                 mutate(seed = x, adjust = "unadjusted", type = "complete"),
#             glm(t2d ~ obese + age + female + cci, data = tmp$amp, family = "binomial") |>
#                 aimTwo::betas_from_mod() |>
#                 filter(predictor == "obese") |>
#                 mutate(seed = x, adjust = "adjusted", type = "complete"),
#             as.data.table(pool1$pooled)[term == "obese", .(predictor = term, beta = estimate)] |>
#                 mutate(seed = x, adjust = "unadjusted", type = "imputed"),
#             as.data.table(pool2$pooled)[term == "obese", .(predictor = term, beta = estimate)] |>
#                 mutate(seed = x, adjust = "adjusted", type = "imputed")
#         )
#     },
#     .progress = TRUE,
#     .options = furrr_options(seed = TRUE)
# )