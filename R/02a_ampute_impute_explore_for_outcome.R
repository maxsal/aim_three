# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice, furrr, glue
)

outcome_phe <- "EM_202.2"

plan(multicore, workers = 10)

# data -------------------------------------------------------------------------
analytic_sub <- qread(glue("data/private/analytic_sub_{outcome_phe}.qs"))
pim <- qread(glue("data/private/mgi/20220822/MGI_PIM0X_20220822_{outcome_phe}.qs"))

# keep these pim vars
pim_sub_vars <- c(
    "id" = "id",
    "case" = "case"
)

cond <- pim |>
    select(any_of(pim_sub_vars))

# merge pim vars into analytic_sub
merged <- analytic_sub |> left_join(cond, by = "id")

# keep limited data
limited <- merged |>
    select(
        age, female, race_eth, cci = charlson_quan_weighted, # covariates
        height = height_cm_med, bmi = bmi_med, smoke = smoke_ever, alcohol = alcohol_ever, # exposures
        case # outcome
    ) |>
    mutate(
        obese = ifelse(bmi >= 30, 1, 0)
    ) |> select(age, female, obese, cci, case)

expanded <- merged |>
    mutate(nhw = as.numeric(race_eth == "NHW")) |>
    select(
        age, female, nhw, cci = charlson_quan_weighted, # covariates
        height = height_cm_med, bmi = bmi_med, smoke = smoke_ever, alcohol = alcohol_ever, # exposures
        case # outcome
    )


expanded$height <- NULL

# MCAR
mcar_ind          <- rbinom(nrow(expanded), 1, 0.5)
expanded$bmi_mcar <- ifelse(mcar_ind == 1, NA, expanded$bmi)

# MAR
make_mar_ind <- function(data, cols, wgts = NULL) {
    if (is.null(wgts)) {
        wgts <- rep(1, length(cols) + 1)
    }
    out <- plogis(
        as.matrix(
            cbind(1, data[, ..cols] |>
                mutate(across(where(is.numeric), ~ scale(.)[, 1])))
    ) %*% wgts) |>
    (\(x) rbinom(length(x), 1, x))()
    message(paste0("Miss%: ", round(mean(out) * 100, 1)))
    return(out)
}

cols <- c("age", "female", "nhw", "cci", "smoke", "alcohol", "case")
mar_prob <- plogis(
    as.matrix(
        cbind(1, expanded[, ..cols] |>
            mutate(across(where(is.numeric), ~ scale(.)[, 1])))
    ) %*% c(-1.5, -1.5, 0.1, -0.1, -0.8, 0.5, 1, 1)
)
mar_ind <- make_mar_ind(
    expanded,
    cols = cols,
    wgts = c(-1.5, -1.5, 0.1, -0.1, -0.8, 0.5, 1)
)
expanded$mar_prob <- mar_prob
expanded$bmi_mar <- ifelse(mar_ind == 1, NA, expanded$bmi)

expanded |>
ggplot(aes(x = mar_prob, color = factor(female))) +
geom_density()

glm(case ~ bmi, data = expanded, family = "binomial")
glm(case ~ bmi_mar, data = expanded, family = "binomial")

# MNAR
cols2 <- c("age", "female", "nhw", "cci", "smoke", "alcohol", "case", "bmi")
mnar_prob <- plogis(
    as.matrix(
        cbind(1, expanded[, ..cols2] |>
            mutate(across(where(is.numeric), ~ scale(.)[, 1])))
    ) %*% c(-1.5, -1.5, 0.1, -0.1, -0.8, 0.5, 1, 1, 1)
)
mnar_ind <- make_mar_ind(
    expanded,
    cols = c(cols, "bmi"),
    wgts = c(-1.5, -1.5, 0.1, -0.1, -0.8, 0.5, 1, 2)
)
expanded$mnar_prob <- mnar_prob
expanded$bmi_mnar <- ifelse(mnar_ind == 1, NA, expanded$bmi)

expanded |>
    ggplot(aes(x = bmi, y = mnar_prob)) +
    geom_point() +
    geom_smooth()

glm(case ~ bmi, data = expanded, family = "binomial")
glm(case ~ bmi_mnar, data = expanded, family = "binomial")

cols3 <- c("age", "female", "nhw", "cci", "smoke", "alcohol", "case", "bmi_mnar")
mnar_imp <- mice(expanded[, ..cols3], m = 10, print = FALSE)
pool(with(mnar_imp, glm(case ~ bmi_mnar, family = "binomial")))


glm(case ~ bmi, data = expanded, family = "binomial") |> summary() |> coef()
glm(case ~ bmi_mcar, data = expanded, family = "binomial") |> summary() |> coef()
glm(case ~ bmi_mar, data = expanded, family = "binomial") |> summary() |> coef()
glm(case ~ bmi_mnar, data = expanded, family = "binomial") |> summary() |> coef()


## impute
# mcar
mcar_cols <- c(cols, "bmi_mcar", "case")
mcar_imp <- mice(expanded[, ..mcar_cols], m = 10, print = FALSE)

# mar
mar_cols <- c(cols, "bmi_mar", "case")
mar_imp <- mice(expanded[, ..mar_cols], m = 10, print = FALSE)

# mnar
mnar_cols <- c(cols, "bmi_mnar", "case")
mnar_imp <- mice(expanded[, ..mnar_cols], m = 10, print = FALSE)

manual_pool <- function(imp, outcome = "case", ex_var = "bmi_mar") {
    l <- rep(NA_real_, imp$m)
    for (i in 1:imp$m) {
        tmp <- complete(imp, i)
        l[i] <- glm(paste0(outcome, " ~ ", ex_var), data = tmp, family = "binomial") |>
            summary() |>
            coef() |>
            (\(y) y[2, 1])()
    }
    return(mean(l, na.rm = TRUE))
}
manual_pool(mar_imp)


summary(pool(with(mcar_imp, glm(case ~ bmi_mcar, family = "binomial"))))
summary(pool(with(mar_imp, glm(case ~ bmi_mar, family = "binomial"))))
summary(pool(with(mnar_imp, glm(case ~ bmi_mnar, family = "binomial"))))

tibble(
    "truth" = glm(case ~ bmi, data = expanded, family = "binomial") |> summary() |> coef() |> (\(y) y[2,1])(),
    "mcar_cc" = glm(case ~ bmi_mcar, data = expanded, family = "binomial") |> summary() |> coef() |> (\(y) y[2,1])(),
    "mcar_imp" = summary(pool(with(mcar_imp, glm(case ~ bmi_mcar, family = "binomial"))))[2, 2],
    "mar_cc" = glm(case ~ bmi_mar, data = expanded, family = "binomial") |> summary() |> coef() |> (\(y) y[2,1])(),
    "mar_imp" = summary(pool(with(mar_imp, glm(case ~ bmi_mar, family = "binomial"))))[2,2],
    "mnar_cc" = glm(case ~ bmi_mnar, data = expanded, family = "binomial") |> summary() |> coef() |> (\(y) y[2,1])(),
    "mnar_imp" = summary(pool(with(mnar_imp, glm(case ~ bmi_mnar, family = "binomial"))))[2,2]
)

mcar_imp <- mice(mcar$amp, m = 10)
mar_imp <- mice(mar$amp, m = 10)
mnar_imp <- mice(mnar$amp, m = 10)

glm(case ~ bmi, data = expanded, family = "binomial") |> summary() |> coef()
glm(case ~ bmi, data = mcar$amp, family = "binomial") |> summary() |> coef()
glm(case ~ bmi, data = mar$amp, family = "binomial") |> summary() |> coef()
glm(case ~ bmi, data = mnar$amp, family = "binomial") |> summary() |> coef()
summary(pool(with(mcar_imp, glm(case ~ bmi, family = "binomial"))))
summary(pool(with(mar_imp, glm(case ~ bmi, family = "binomial"))))
summary(pool(with(mnar_imp, glm(case ~ bmi, family = "binomial"))))

test_manual_missing <- map_dfr(
    1:25,
    \(x) {
        set.seed(x)

        # MCAR
        mcar_ind <- rbinom(nrow(expanded), 1, 0.3)
        expanded$bmi_mcar <- ifelse(mcar_ind == 1, NA, expanded$bmi)

        cols <- c("age", "female", "nhw", "cci", "smoke", "alcohol")
        # MAR

        # try_these_weights <- function(wgts) {
        #     tmp_probs <- plogis(
        #         as.matrix(
        #             cbind(1, expanded[, ..cols] |>
        #                 mutate(across(where(is.numeric), ~ scale(.)[, 1])))
        #         ) %*% wgts
        #     )
        #     tmp_ind <- rbinom(n = length(tmp_probs), size = 1, prob = tmp_probs)

        #     tibble(
        #         x = tmp_probs
        #     ) |>
        #         ggplot(aes(x = x)) +
        #         geom_histogram(bins = 30) +
        #         geom_text(
        #             aes(x = 0.5, y = Inf),
        #             hjust = 0.5, vjust = 2, label = round(mean(tmp_ind) * 100, 1)
        #         ) +
        #         ms::theme_ms()
        # }
        # try_these_weights(
        #     c(-1.5, -1.5, 0.1, -0.1, -0.8, 0.5, 1)
        # )

        mar_ind <- plogis(
            as.matrix(
                cbind(1, expanded[, ..cols] |>
                    mutate(across(where(is.numeric), ~ scale(.)[, 1])))
            ) %*% c(-1.5, -1, 0.1, -0.1, -0.8, 0.5, 1)
        ) |>
        (\(x) rbinom(n = length(x), size = 1, prob = x))()
        expanded$bmi_mar <- ifelse(mar_ind == 1, NA, expanded$bmi)

        # MNAR
        cols_bmi <- c(cols, "bmi")
        mnar_ind <- plogis(
            as.matrix(
                cbind(1, expanded[, ..cols_bmi] |>
                    mutate(across(where(is.numeric), ~ scale(.)[, 1])))
            ) %*% c(0, 1, -1, 2, 1, -2, 1, 2)
        ) |> (\(x) rbinom(n = length(x), size = 1, prob = x))()
        expanded$bmi_mnar <- ifelse(mnar_ind == 1, NA, expanded$bmi)

        truth     <- glm(case ~ bmi, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi") |> pull(beta)
        mcar_beta <- glm(case ~ bmi_mcar, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi_mcar") |> pull(beta)
        mar_beta  <- glm(case ~ bmi_mar, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi_mar") |> pull(beta)
        mnar_beta <- glm(case ~ bmi_mnar, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi_mnar") |> pull(beta)

        # impute
        ## MCAR
        mcar_vars <- c(cols, "bmi_mcar")
        case <- expanded$case
        mcar_imp <- mice(expanded[, ..mcar_vars], m = 10, print = FALSE)
        mcar_imp_beta <- rep(NA_real_, 10)
        for (i in 1:10) {
            tmp <- complete(mcar_imp, i)
            mcar_imp_beta[i] <- glm(case ~ bmi_mcar, data = tmp, family = "binomial") |>
                summary() |>
                coef() |>
                (\(y) y[2, 1])()
        }
        mcar_imp_beta <- mean(mcar_imp_beta, na.rm = TRUE)

        ## MAR
        mar_vars <- c(cols, "bmi_mar")
        mar_imp <- mice(expanded[, ..mar_vars], m = 10, print = FALSE)
        # mar_imp$data |>
        #     ggplot(aes(x = bmi_mar)) +
        #     geom_density() +
        #     geom_density(data = complete(mar_imp, 1), aes(x = bmi_mar), color = "red")

        # pool(with(mar_imp, glm(case ~ bmi_mar, family = "binomial")))


        # glm(case ~ bmi, data = expanded, family = "binomial")
        # glm(case ~ bmi_mar, data = complete(mar_imp, 1), family = "binomial")

        mar_imp_beta <- rep(NA_real_, 10)
        for (i in 1:10) {
            tmp <- complete(mar_imp, i)
            mar_imp_beta[i] <- glm(case ~ bmi_mar, data = tmp, family = "binomial") |>
                summary() |>
                coef() |>
                (\(y) y[2, 1])()
        }
        mar_imp_beta <- mean(mar_imp_beta, na.rm = TRUE)
        
        ## MNAR
        mnar_vars <- c(cols, "bmi_mnar")
        mnar_imp <- mice(expanded[, ..mnar_vars], m = 10, print = FALSE)
        # mnar_imp_beta <- pool(with(mnar_imp, glm(case ~ bmi_mnar, family = "binomial"))) |>
        #     summary() |>
        #     as_tibble() |>
        #     filter(term == "bmi_mnar") |>
        #     pull("estimate")

        mnar_imp_beta <- rep(NA_real_, 10)
        for (i in 1:10) {
            tmp <- complete(mnar_imp, i)
            mnar_imp_beta[i] <- glm(case ~ bmi_mnar, data = tmp, family = "binomial") |>
                summary() |>
                coef() |>
                (\(y) y[2, 1])()
        }
        mnar_imp_beta <- mean(mnar_imp_beta, na.rm = TRUE)

        tibble(
            mech = c("mcar", "mar", "mnar", "mcar", "mar", "mnar"),
            type = c("cc", "cc", "cc", "imp", "imp", "imp"),
            bias = c(
                mcar_beta - truth,
                mar_beta - truth,
                mnar_beta - truth,
                mcar_imp_beta - truth,
                mar_imp_beta - truth,
                mnar_imp_beta - truth
            ),
            seed = x
        )
    },
    .progress = TRUE
)

test_manual_missing |>
    ggplot(aes(x = bias, color = mech, linetype = type)) +
    geom_vline(xintercept = 0) +
    geom_density(linewidth = 1) +
    facet_wrap(~type, ncol = 1, scales = "free") +
    ms::theme_ms()


# test_mar <- cascade_amp_imp(
#     data = expanded,
#     outcome = "case",
#     predictor = "bmi",
#     covs = c("age", "female", "nhw", "cci"),
#     exposures = c("smoke", "alcohol"),
#     patterns = NULL,
#     freq  = NULL,
#     weights = NULL,
#     miss_prop = 0.3,
#     mech = "MAR",
#     n_seeds = 100,
#     n_imp = 10,
#     .progress = TRUE
# )

# truths <- tibble(
#     adjust = c("un_adj", "cov_adj", "full_adj"),
#     beta = c(
#         glm(case ~ bmi, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi") |> pull(beta),
#         glm(case ~ bmi + age + female + nhw + cci, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi") |> pull(beta),
#         glm(case ~ bmi + age + female + nhw + cci + smoke + alcohol, data = expanded, family = "binomial") |> aimTwo::betas_from_mod() |> filter(predictor == "bmi") |> pull(beta)
#     )
# )

# test_mar |>
#     ggplot(aes(x = beta, linetype = type, color = type)) +
#     geom_vline(data = truths, aes(xintercept = beta)) +
#     geom_density(linewidth = 1) +
#     facet_wrap(~adjust, ncol = 1)


# l <- rep(0, 10)
# for (i in 1:10) {
#     x = complete(mar_imp, i)
#     l[i] <- glm(case ~ bmi_mar, data = x, family = "binomial") |> summary() |> coef() |> (\(y) y[2,1])()
# }

# explore ampute base patterns
base_pattern <- ampute(expanded)$patterns
base_pattern$age <- base_pattern$female <- base_pattern$nhw <-  base_pattern$cci <- base_pattern$case <- 1
keep <- apply(base_pattern, 1, \(x) length(unique(x[!is.na(x)])) != 1)
base_pattern[keep, ]

new_patterns <- rbind(
    base_pattern[keep, ],
    c(1, 1, 1, 1, 0, 0, 1, 1),
    c(1, 1, 1, 1, 1, 0, 0, 1),
    c(1, 1, 1, 1, 0, 1, 0, 1),
    c(1, 1, 1, 1, 0, 0, 0, 1)
)

result <- ampute(expanded, patterns = new_patterns)

patt <- rep(1, ncol(limited))
patt[which(names(limited) == "obese")] <- 0
wts <- rep(1, ncol(limited))
wts[which(names(limited) %in% c("obese", "case"))] <- 0

# map_dbl(mar_0.5$amp, \(x) sum(is.na(x))) |> enframe() |> arrange(desc(value))

# setup for ampute/impute function
outcome   <- "case"
predictor <- .predictor <- "bmi"
covs      <- c("age", "female", "nhw", "cci")
exposures <- c("bmi", "smoke", "alcohol")
exposures <- exposures[!exposures %in% predictor]
n_seeds   <- 100
n_imps    <- 10
miss_prop <- 0.3

# complete case
# unadjusted
un_adj_cc <- glm(paste0(outcome, "~", predictor), data = expanded, family = "binomial") |> aimTwo::betas_from_mod()
# covariate adjusted
cov_adj_cc <- glm(paste0(outcome, "~", paste0(c(predictor, covs), collapse = "+")), data = expanded, family = "binomial") |> aimTwo::betas_from_mod()
# full adjusted
full_adj_cc <- glm(paste0(outcome, "~", paste0(c(predictor, covs, exposures), collapse = "+")), data = expanded, family = "binomial") |> aimTwo::betas_from_mod()

vlines <- tibble(
    adjust = c("Unadjusted", "Covariate adjusted", "Fully adjusted"),
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
    outcome = outcome,
    predictor = predictor,
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
    outcome = outcome,
    predictor = predictor,
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
    outcome   = outcome,
    predictor = predictor,
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
    mcar |> mutate(`Missingness Mechanism` = "MCAR"),
    mar |> mutate(`Missingness Mechanism` = "MAR"),
    mnar |> mutate(`Missingness Mechanism` = "MNAR")
) |>
    mutate(
    adjust = factor(case_when(
        adjust == "un_adj" ~ "Unadjusted",
        adjust == "cov_adj" ~ "Covariate adjusted",
        adjust == "full_adj" ~ "Fully adjusted"
    ), levels = c("Unadjusted", "Covariate adjusted", "Fully adjusted")),
    `Missingness Mechanism` = factor(
        `Missingness Mechanism`,
        levels = c("MCAR", "MAR", "MNAR"))
    ) |>
    ggplot(aes(x = beta, color = `Missingness Mechanism`, linetype = type)) +
    geom_vline(data = vlines, aes(xintercept = beta)) +
    # geom_histogram(bins = 30, alpha = 0.5, position = "dodge") +
    geom_density(linewidth = 1, alpha = 0.5) +
    facet_wrap(~adjust, ncol = 1) +
    scale_color_ms() +
    theme_ms()

ggsave(
    glue("bin/test_plot_{outcome_phe}.pdf"),
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