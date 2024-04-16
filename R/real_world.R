ms::libri(ms, qs, tidyverse, data.table, survey, mice, cli, glue, patchwork)

# demo
d <- qread("data/private/analytic_full.qs")
# pim
pim <- qread("data/private/mgi_pim0x_20230322.qs")
# prs
ex_prs <- qread("data/private/ex_prs_processed.qs")
prs_sub_vars <- c(
    "id" = "id",
    "bmi_prs" = "bmi",
    "t2d_prs" = "t2d",
    "glucose_prs" = "glucose"
)
prs <- ex_prs |>
    mutate(
        bmi_prs = scale(bmi)[, 1],
        t2d_prs = scale(t2d)[, 1],
        glucose_prs = scale(glucose)[, 1]
    ) |>
    select(id, bmi_prs, t2d_prs, glucose_prs)

# weights
w <- qread("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/mgi/20220822/weightsx_20220822_comb.qs")

# subset
d   <- d[id %in% w[!is.na(ip_selection), id], ]
pim <- pim[id %in% unique(d[, id]), ]
w   <- w[id %in% unique(d[, id]), .(id, weights = ip_selection)]
prs <- prs[id %in% d[, id], ]

d_sub <- d |>
    left_join(
        pim |>
            select(id, t2d = EM_202.2),
        by = "id"
    ) |>
    mutate(nhw = case_when(
        race_eth == "NHW" ~ 1,
        race_eth != "NHW" ~ 0,
        TRUE ~ NA_real_
    )) |>
    select(id, age, female, nhw, bmi = bmi_med, smoke = smoke_ever, t2d, glucose = glucose_med)    
d_sub_id <- d_sub |> select(id)
d_sub <- d_sub |> select(-id)

# missingness plot
map_dfr(names(d),
    \(x) {
        data.table(
            var = x,
            miss_prop = mean(is.na(d[[x]]))
        )
    }) |>
    arrange(desc(miss_prop)) |>
    filter(miss_prop > 0) |>
    ggplot(aes(x = reorder(var, miss_prop), y = miss_prop)) +
    geom_segment(aes(xend = var, yend = 0), linewidth = 0.5) +
    geom_point(size = 3) +
    coord_flip() +
    labs(title = "Missingness by Variable",
         x = "",
         y = "Proportion Missing") +
    scale_y_continuous(labels = scales::percent) +
    ms::theme_ms() +
    theme(
        axis.text.y = element_text(size = 6),
        panel.grid.major.y = element_blank()
    )

# functions --------------------------------------------------------------------
# return nice results from glm estimates
.glm <- function(f, data, var, fam = "gaussian") {
    mod <- glm(f, data = data, family = fam)
    est <- mod$coefficients[[var]]
    ci <- suppressMessages({
        confint(mod)
    })[var, ]
    data.table(
        var = var,
        est = est,
        lo = ci[[1]],
        hi = ci[[2]]
    )
}

# return nice results from svyglm estimates
.svyglm <- function(f, design, var, fam = "gaussian") {
    mod <- svyglm(as.formula(f), design = design, family = fam)
    est <- coef(mod)[[var]]
    ci <- confint(mod)[var, ]
    variance <- vcov(mod)[var, var]
    data.table(
        variable = var,
        est      = est,
        lo       = ci[[1]],
        hi       = ci[[2]],
        var      = variance
    )
}
# return nice pooled results from imputed glm estimates
.glm_pool <- function(f, imp, var, fam = "gaussian") {
    mod <- map(
        seq_len(imp$m), 
        \(i) {
            glm(as.formula(f), data = complete(imp, i), family = fam)
        }
    ) |>
    pool() |>
    summary(conf.int = TRUE)

    row_ind <- which(mod$term == var)
    col_ind <- which(names(mod) %in% c("estimate", "2.5 %", "97.5 %"))
    data.table(
        var = var,
        est = mod[row_ind, col_ind[1]],
        lo = mod[row_ind, col_ind[2]],
        hi = mod[row_ind, col_ind[3]]
    )
}

.glm_pool_obese1 <- function(f, imp, var, fam = "binomial") {
    mod <- map(
        seq_len(imp$m),
        \(i) {
            glm(
                as.formula(f),
                data = complete(imp, i) |>
                    mutate(
                        obese1 = case_when(
                            bmi >= 30 & bmi < 35 ~ 1,
                            bmi >= 18.5 & bmi < 25 ~ 0,
                            TRUE ~ NA_real_
                        )
                    ) |>
                    drop_na(obese1) |>
                    select(-bmi),
                family = fam)
        }
    ) |>
        pool() |>
        summary(conf.int = TRUE)

    row_ind <- which(mod$term == var)
    col_ind <- which(names(mod) %in% c("estimate", "2.5 %", "97.5 %"))
    data.table(
        var = var,
        est = mod[row_ind, col_ind[1]],
        lo = mod[row_ind, col_ind[2]],
        hi = mod[row_ind, col_ind[3]]
    )
}

.glm_pool_obese2 <- function(f, imp, var, fam = "binomial") {
    mod <- map(
        seq_len(imp$m),
        \(i) {
            glm(
                as.formula(f),
                data = complete(imp, i) |>
                    mutate(
                        obese2 = case_when(
                            bmi >= 35 & bmi < 40 ~ 1,
                            bmi >= 18.5 & bmi < 25 ~ 0,
                            TRUE ~ NA_real_
                        )
                    ) |>
                    drop_na(obese2) |>
                    select(-bmi),
                family = fam
            )
        }
    ) |>
        pool() |>
        summary(conf.int = TRUE)

    row_ind <- which(mod$term == var)
    col_ind <- which(names(mod) %in% c("estimate", "2.5 %", "97.5 %"))
    data.table(
        var = var,
        est = mod[row_ind, col_ind[1]],
        lo = mod[row_ind, col_ind[2]],
        hi = mod[row_ind, col_ind[3]]
    )
}

# pooling survey-weighted glm estimates from a MICE object
.svyglm_pool <- function(f, imp_object, weights, var, fam = "gaussian") {
    svyglm_res <- map_dfr(
        seq_len(imp_object$m),
        \(i) {
            impw <- cbind(complete(imp_object, i), weights = weights)
            impw_dsn <- svydesign(
                ids     = ~1,
                data    = impw,
                weights = ~weights
            )
            .svyglm(f, var = var, design = impw_dsn, fam = fam)
        }
    )

    pooled_mean <- mean(svyglm_res$est)
    within_var <- mean(svyglm_res$var)
    between_var <- var(svyglm_res$est)
    df_correction <- (nrow(svyglm_res) + 1) / nrow(svyglm_res)
    pooled_var <- within_var + (between_var * df_correction)
    pooled_se <- sqrt(pooled_var)
    lo <- pooled_mean - (qnorm(0.975) * pooled_se)
    hi <- pooled_mean + (qnorm(0.975) * pooled_se)

    data.table(var = var, est = pooled_mean, lo = lo, hi = hi, var = pooled_var)
}

.svyglm_pool_obese1 <- function(f, imp_object, weights, var, fam = "gaussian") {
    svyglm_res <- map_dfr(
        seq_len(imp_object$m),
        \(i) {
            impw <- cbind(complete(imp_object, i), weights = weights) |>
                mutate(
                    obese1 = case_when(
                        bmi >= 30 & bmi < 35 ~ 1,
                        bmi >= 18.5 & bmi < 25 ~ 0,
                        TRUE ~ NA_real_
                    )
                ) |>
                drop_na(obese1) |>
                select(-bmi)
            impw_dsn <- svydesign(
                ids     = ~1,
                data    = impw,
                weights = ~weights
            )
            .svyglm(f, var = var, design = impw_dsn, fam = fam)
        }
    )

    pooled_mean <- mean(svyglm_res$est)
    within_var <- mean(svyglm_res$var)
    between_var <- var(svyglm_res$est)
    df_correction <- (nrow(svyglm_res) + 1) / nrow(svyglm_res)
    pooled_var <- within_var + (between_var * df_correction)
    pooled_se <- sqrt(pooled_var)
    lo <- pooled_mean - (qnorm(0.975) * pooled_se)
    hi <- pooled_mean + (qnorm(0.975) * pooled_se)

    data.table(var = var, est = pooled_mean, lo = lo, hi = hi, var = pooled_var)
}

.svyglm_pool_obese2 <- function(f, imp_object, weights, var, fam = "gaussian") {
    svyglm_res <- map_dfr(
        seq_len(imp_object$m),
        \(i) {
            impw <- cbind(complete(imp_object, i), weights = weights) |>
                mutate(
                    obese2 = case_when(
                        bmi >= 35 & bmi < 40 ~ 1,
                        bmi >= 18.5 & bmi < 25 ~ 0,
                        TRUE ~ NA_real_
                    )
                ) |>
                drop_na(obese2) |>
                select(-bmi)
            impw_dsn <- svydesign(
                ids     = ~1,
                data    = impw,
                weights = ~weights
            )
            .svyglm(f, var = var, design = impw_dsn, fam = fam)
        }
    )

    pooled_mean <- mean(svyglm_res$est)
    within_var <- mean(svyglm_res$var)
    between_var <- var(svyglm_res$est)
    df_correction <- (nrow(svyglm_res) + 1) / nrow(svyglm_res)
    pooled_var <- within_var + (between_var * df_correction)
    pooled_se <- sqrt(pooled_var)
    lo <- pooled_mean - (qnorm(0.975) * pooled_se)
    hi <- pooled_mean + (qnorm(0.975) * pooled_se)

    data.table(var = var, est = pooled_mean, lo = lo, hi = hi, var = pooled_var)
}

# GLUCOSE ~ BMI ----------------------------------------------------------------
outcome     <- "glucose"
exposure    <- "bmi"
covariates  <- c("age", "female", "nhw", "smoke")
d_analytic  <- d_sub |>
    select(all_of(c(outcome, exposure, covariates)))
non_outcome <- names(d_analytic) |> setdiff(outcome)
f           <- paste0(outcome, " ~ ", exposure)
f_cov       <- paste0(outcome, " ~ ", paste(non_outcome, collapse = " + "))

# COMPLETE CASE
# standard
cc_un_est  <- .glm(f, data = d_analytic, var = exposure)
cc_adj_est <- .glm(f_cov, data = d_analytic, var = exposure)

# with weights
d_analytic_w <- cbind(d_sub_id, d_analytic) |>
    left_join(w, by = "id")
cc_dsn <- svydesign(
    ids = ~1,
    data = d_analytic_w,
    weights = ~weights
)

ccw_un_est  <- .svyglm(f, var = exposure, design = cc_dsn)
ccw_adj_est <- .svyglm(f_cov, var = exposure, design = cc_dsn)

# IMPUTATION
imp     <- mice::mice(d_analytic, m = 5, seed = 123, printFlag = FALSE)

imp_un_est  <- .glm_pool(f = f, imp = imp, var = exposure)
imp_adj_est <- .glm_pool(f = f_cov, imp = imp, var = exposure)

impw_un_est  <- .svyglm_pool(f, imp, weights = w[, weights], var = exposure)
impw_adj_est <- .svyglm_pool(f_cov, imp, weights = w[, weights], var = exposure)

# IMPUTATION w PRS
imp_prs <- mice::mice(cbind(d_analytic, prs[, .(glucose_prs, bmi_prs)]), m = 5, seed = 123, printFlag = FALSE)

imp_prs_un_est <- .glm_pool(f = f, imp = imp_prs, var = exposure)
imp_prs_adj_est <- .glm_pool(f = f_cov, imp = imp_prs, var = exposure)

impw_prs_un_est <- .svyglm_pool(f, imp_prs, weights = w[, weights], var = exposure)
impw_prs_adj_est <- .svyglm_pool(f_cov, imp_prs, weights = w[, weights], var = exposure)

glu_bmi_p_nhanes_est <- tibble(
    adjust = factor(c("Unadjusted", "Adjusted"), c("Unadjusted", "Adjusted")),
    xmin = c(-Inf, -Inf),
    xmax = c(Inf, Inf),
    ymin = c(0.624, 0.601),
    ymax = c(1.098, 1.061),
    est = c(NA_real_, NA_real_),
    lo = c(NA_real_, NA_real_),
    hi = c(NA_real_, NA_real_),
    mdm = c(NA_character_, NA_character_),
    weight = c(NA_character_, NA_character_)
)

gluc_bmi_plot <- tibble(
    mdm    = c(rep("Complete case", 4), rep ("Imputed", 4), rep("PRS-informed", 4)), # missing data method
    weight = c(rep(c("Unweighted", "Unweighted", "Weighted", "Weighted"), 3)),
    adjust = factor(c(rep(c("Unadjusted", "Adjusted"), 6)), c("Unadjusted", "Adjusted")),
    est    = c(
        cc_un_est$est, cc_adj_est$est, ccw_un_est$est, ccw_adj_est$est,
        imp_un_est$est, imp_adj_est$est, impw_un_est$est, impw_adj_est$est,
        imp_prs_un_est$est, imp_prs_adj_est$est, impw_prs_un_est$est, impw_prs_adj_est$est
    ),
    lo    = c(
        cc_un_est$lo, cc_adj_est$lo, ccw_un_est$lo, ccw_adj_est$lo,
        imp_un_est$lo, imp_adj_est$lo, impw_un_est$lo, impw_adj_est$lo,
        imp_prs_un_est$lo, imp_prs_adj_est$lo, impw_prs_un_est$lo, impw_prs_adj_est$lo
    ),
    hi    = c(
        cc_un_est$hi, cc_adj_est$hi, ccw_un_est$hi, ccw_adj_est$hi,
        imp_un_est$hi, imp_adj_est$hi, impw_un_est$hi, impw_adj_est$hi,
        imp_prs_un_est$hi, imp_prs_adj_est$hi, impw_prs_un_est$hi, impw_prs_adj_est$hi
    )
) |>
    ggplot(aes(x = mdm, y = est, color = weight, shape = weight)) +
    geom_rect(
        data = glu_bmi_p_nhanes_est,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = NULL, y = NULL),
        fill = "grey",
        alpha = 0.2,
        show.legend = FALSE
        # inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0) +
    geom_pointrange(
        aes(ymin = lo, ymax = hi),
        position = position_dodge(width = 0.2),
        size = 1,
        linewidth = 1) +
    # geom_point(position = position_dodge(width = 0.2), size = 3) +
    coord_flip() +
    facet_wrap(~adjust) +
    labs(
        title = "BMI coefficient for glucose by adjustment and missing data method",
        x = "Missing data method",
        y = "Beta coefficient (95% CI)",
        caption = "Shaded region represents estimated range using NHANES 2017-2020 (P) data (fasting glucose)"
    ) +
    ms::scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.title = element_blank()
    )



# T2D ~ BMI --------------------------------------------------------------------
outcome    <- "t2d"
family     <- "binomial"
exposure   <- "bmi"
covariates <- c("age", "female", "nhw", "smoke")
d_analytic <- d_sub |>
    select(all_of(c(outcome, exposure, covariates)))
non_outcome <- names(d_analytic) |> setdiff(outcome)
f <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(non_outcome, collapse = " + "))

# COMPLETE CASE
# standard
cc_un_est <- .glm(f, data = d_analytic, var = exposure, fam = family)
cc_adj_est <- .glm(f_cov, data = d_analytic, var = exposure, fam = family)

# with weights
d_analytic_w <- cbind(d_sub_id, d_analytic) |>
    left_join(w, by = "id")
cc_dsn <- svydesign(
    ids = ~1,
    data = d_analytic_w,
    weights = ~weights
)

svyglm_family <- switch(
    family,
    "binomial" = "quasibinomial",
    "guassian" = "gaussian"
)
ccw_un_est <- .svyglm(f, var = exposure, design = cc_dsn, fam = svyglm_family)
ccw_adj_est <- .svyglm(f_cov, var = exposure, design = cc_dsn, fam = svyglm_family)

# IMPUTATION
imp <- mice::mice(d_analytic, m = 5, seed = 123, printFlag = FALSE)

imp_un_est  <- .glm_pool(f = f, imp = imp, var = exposure, fam = family)
imp_adj_est <- .glm_pool(f = f_cov, imp = imp, var = exposure, fam = family)

impw_un_est <- .svyglm_pool(f, imp, weights = w[, weights], var = exposure, fam = svyglm_family)
impw_adj_est <- .svyglm_pool(f_cov, imp, weights = w[, weights], var = exposure, fam = svyglm_family)

# IMPUTATION w PRS
imp_prs <- mice::mice(cbind(d_analytic, prs[, .(bmi_prs, t2d_prs)]), m = 5, seed = 123, printFlag = FALSE)

imp_prs_un_est <- .glm_pool(f = f, imp = imp_prs, var = exposure)
imp_prs_adj_est <- .glm_pool(f = f_cov, imp = imp_prs, var = exposure)

impw_prs_un_est <- .svyglm_pool(f, imp_prs, weights = w[, weights], var = exposure)
impw_prs_adj_est <- .svyglm_pool(f_cov, imp_prs, weights = w[, weights], var = exposure)

t2d_bmi_plot <- tibble(
    mdm    = c(rep("Complete case", 4), rep ("Imputed", 4), rep("PRS-informed", 4)), # missing data method
    weight = c(rep(c("Unweighted", "Unweighted", "Weighted", "Weighted"), 3)),
    adjust = c(rep(c("Unadjusted", "Adjusted"), 6)),
    est    = exp(c(
        cc_un_est$est, cc_adj_est$est, ccw_un_est$est, ccw_adj_est$est,
        imp_un_est$est, imp_adj_est$est, impw_un_est$est, impw_adj_est$est,
        imp_prs_un_est$est, imp_prs_adj_est$est, impw_prs_un_est$est, impw_prs_adj_est$est
    )),
    lo    = exp(c(
        cc_un_est$lo, cc_adj_est$lo, ccw_un_est$lo, ccw_adj_est$lo,
        imp_un_est$lo, imp_adj_est$lo, impw_un_est$lo, impw_adj_est$lo,
        imp_prs_un_est$lo, imp_prs_adj_est$lo, impw_prs_un_est$lo, impw_prs_adj_est$lo
    )),
    hi    = exp(c(
        cc_un_est$hi, cc_adj_est$hi, ccw_un_est$hi, ccw_adj_est$hi,
        imp_un_est$hi, imp_adj_est$hi, impw_un_est$hi, impw_adj_est$hi,
        imp_prs_un_est$hi, imp_prs_adj_est$hi, impw_prs_un_est$hi, impw_prs_adj_est$hi
    ))
) |>
    ggplot(aes(x = mdm, y = est, color = weight, shape = weight)) +
    geom_rect(
        aes(xmin = -Inf, xmax = Inf, ymin = 1.1, ymax = 1.13),
        fill = "grey",
        alpha = 0.1,
        inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 1) +
    geom_pointrange(
        aes(ymin = lo, ymax = hi),
        position = position_dodge(width = 0.2),
        size = 1,
        linewidth = 1
    ) +
    # geom_point(position = position_dodge(width = 0.2), size = 3) +
    coord_flip() +
    facet_wrap(~ factor(adjust, c("Unadjusted", "Adjusted"))) +
    labs(
        title = "BMI coefficient for type 2 diabetes by adjustment and missing data method",
        x = "Missing data method",
        y = "Odds ratio (95% CI)",
        caption = "Shaded region represents (covariate-adjusted) estimated range reported by Strings et al. (2023)"
    ) +
    ms::scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.title = element_blank()
    )

# OBESITY (class 1: [30,35)]) ~ T2D --------------------------------------------
outcome <- "obese1"
family  <- "binomial"
exposure <- "t2d"
covariates <- c("age", "female", "nhw", "smoke")
d_analytic <- cbind(d_sub_id, d_sub) |>
    mutate(obese1 = case_when(
        bmi >= 30 & bmi < 35 ~ 1,
        bmi >= 18.5 & bmi < 25 ~ 0,
        TRUE ~ NA_real_
    )) |>
    drop_na(obese1) |>
    select(all_of(c("id", outcome, exposure, covariates)))
d_analytic_id <- d_analytic |> select(id)
d_analytic <- d_analytic |> select(-id)

# retain second dataset with raw BMI to impute prior to generating obesity variable
d_analytic2 <- d_sub |>
    select(all_of(c("bmi", exposure, covariates)))
non_outcome <- names(d_analytic) |> setdiff(outcome)
f <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(non_outcome, collapse = " + "))

# COMPLETE CASE
# standard
cc_un_est  <- .glm(f, data = d_analytic, var = exposure, fam = family)
cc_adj_est <- .glm(f_cov, data = d_analytic, var = exposure, fam = family)

# with weights
d_analytic_w <- cbind(d_analytic_id, d_analytic) |>
    left_join(w, by = "id")
cc_dsn <- svydesign(
    ids = ~1,
    data = d_analytic_w,
    weights = ~weights
)

svyglm_family <- switch(family,
    "binomial" = "quasibinomial",
    "guassian" = "gaussian"
)
ccw_un_est  <- .svyglm(f, var = exposure, design = cc_dsn, fam = svyglm_family)
ccw_adj_est <- .svyglm(f_cov, var = exposure, design = cc_dsn, fam = svyglm_family)

# IMPUTATION
imp <- mice::mice(d_analytic2, m = 5, seed = 123, printFlag = FALSE) # impute BMI before generating variable

imp_un_est  <- .glm_pool_obese1(f = f, imp = imp, var = exposure, fam = family)
imp_adj_est <- .glm_pool_obese1(f = f_cov, imp = imp, var = exposure, fam = family)

impw_un_est  <- .svyglm_pool_obese1(f, imp, weights = w[, weights], var = exposure, fam = svyglm_family)
impw_adj_est <- .svyglm_pool_obese1(f_cov, imp, weights = w[, weights], var = exposure, fam = svyglm_family)

# IMPUTATION + PRS
imp_prs <- mice::mice(cbind(d_analytic2, prs[, .(bmi_prs, t2d_prs)]), m = 5, seed = 123, printFlag = FALSE) # impute BMI before generating variable

imp_prs_un_est <- .glm_pool_obese1(f = f, imp = imp_prs, var = exposure, fam = family)
imp_prs_adj_est <- .glm_pool_obese1(f = f_cov, imp = imp_prs, var = exposure, fam = family)

impw_prs_un_est <- .svyglm_pool_obese1(f, imp_prs, weights = w[, weights], var = exposure, fam = svyglm_family)
impw_prs_adj_est <- .svyglm_pool_obese1(f_cov, imp_prs, weights = w[, weights], var = exposure, fam = svyglm_family)


ob1_t2d_plot <- tibble(
    mdm    = c(rep("Complete case", 4), rep ("Imputed", 4), rep("PRS-informed", 4)), # missing data method
    weight = c(rep(c("Unweighted", "Unweighted", "Weighted", "Weighted"), 3)),
    adjust = c(rep(c("Unadjusted", "Adjusted"), 6)),
    est    = exp(c(
        cc_un_est$est, cc_adj_est$est, ccw_un_est$est, ccw_adj_est$est,
        imp_un_est$est, imp_adj_est$est, impw_un_est$est, impw_adj_est$est,
        imp_prs_un_est$est, imp_prs_adj_est$est, impw_prs_un_est$est, impw_prs_adj_est$est
    )),
    lo    = exp(c(
        cc_un_est$lo, cc_adj_est$lo, ccw_un_est$lo, ccw_adj_est$lo,
        imp_un_est$lo, imp_adj_est$lo, impw_un_est$lo, impw_adj_est$lo,
        imp_prs_un_est$lo, imp_prs_adj_est$lo, impw_prs_un_est$lo, impw_prs_adj_est$lo
    )),
    hi    = exp(c(
        cc_un_est$hi, cc_adj_est$hi, ccw_un_est$hi, ccw_adj_est$hi,
        imp_un_est$hi, imp_adj_est$hi, impw_un_est$hi, impw_adj_est$hi,
        imp_prs_un_est$hi, imp_prs_adj_est$hi, impw_prs_un_est$hi, impw_prs_adj_est$hi
    ))
) |>
    ggplot(aes(x = mdm, y = est, color = weight, shape = weight)) +
    geom_rect(
        aes(xmin = -Inf, xmax = Inf, ymin = 2.3, ymax = 2.6),
        fill = "grey",
        alpha = 0.1,
        inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 1) +
    geom_pointrange(
        aes(ymin = lo, ymax = hi),
        position = position_dodge(width = 0.2),
        size = 1,
        linewidth = 1
    ) +
    # geom_point(position = position_dodge(width = 0.2), size = 3) +
    coord_flip() +
    facet_wrap(~ factor(adjust, c("Unadjusted", "Adjusted"))) +
    labs(
        title = "Type 2 diabetes coefficient for class-1 obesity by adjustment and missing data method",
        x = "Missing data method",
        y = "Odds ratio (95% CI)",
        caption = "Shaded region represents (covariate-adjusted) estimated range reported by Ganz et al. (2014)"
    ) +
    ms::scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.title = element_blank()
    )

# OBESITY (class 2: [35,40)]) ~ T2D --------------------------------------------
outcome <- "obese2"
family <- "binomial"
exposure <- "t2d"
covariates <- c("age", "female", "nhw", "smoke")
d_analytic <- cbind(d_sub_id, d_sub) |>
    mutate(obese2 = case_when(
        bmi >= 35 & bmi < 40 ~ 1,
        bmi >= 18.5 & bmi < 25 ~ 0,
        TRUE ~ NA_real_
    )) |>
    drop_na(obese2) |>
    select(all_of(c("id", outcome, exposure, covariates)))
d_analytic_id <- d_analytic |> select(id)
d_analytic <- d_analytic |> select(-id)

# retain second dataset with raw BMI to impute prior to generating obesity variable
d_analytic2 <- d_sub |>
    select(all_of(c("bmi", exposure, covariates)))
non_outcome <- names(d_analytic) |> setdiff(outcome)
f <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(non_outcome, collapse = " + "))

# COMPLETE CASE
# standard
cc_un_est <- .glm(f, data = d_analytic, var = exposure, fam = family)
cc_adj_est <- .glm(f_cov, data = d_analytic, var = exposure, fam = family)

# with weights
d_analytic_w <- cbind(d_analytic_id, d_analytic) |>
    left_join(w, by = "id")
cc_dsn <- svydesign(
    ids = ~1,
    data = d_analytic_w,
    weights = ~weights
)

svyglm_family <- switch(family,
    "binomial" = "quasibinomial",
    "guassian" = "gaussian"
)
ccw_un_est <- .svyglm(f, var = exposure, design = cc_dsn, fam = svyglm_family)
ccw_adj_est <- .svyglm(f_cov, var = exposure, design = cc_dsn, fam = svyglm_family)

# IMPUTATION
imp <- mice::mice(d_analytic2, m = 5, seed = 123, printFlag = FALSE) # impute BMI before generating variable

imp_un_est <- .glm_pool_obese2(f = f, imp = imp, var = exposure, fam = family)
imp_adj_est <- .glm_pool_obese2(f = f_cov, imp = imp, var = exposure, fam = family)

impw_un_est <- .svyglm_pool_obese2(f, imp, weights = w[, weights], var = exposure, fam = svyglm_family)
impw_adj_est <- .svyglm_pool_obese2(f_cov, imp, weights = w[, weights], var = exposure, fam = svyglm_family)

# IMPUTATION + PRS
imp_prs <- mice::mice(cbind(d_analytic2, prs[, .(bmi_prs, t2d_prs)]), m = 5, seed = 123, printFlag = FALSE) # impute BMI before generating variable

imp_prs_un_est <- .glm_pool_obese2(f = f, imp = imp_prs, var = exposure, fam = family)
imp_prs_adj_est <- .glm_pool_obese2(f = f_cov, imp = imp_prs, var = exposure, fam = family)

impw_prs_un_est <- .svyglm_pool_obese2(f, imp_prs, weights = w[, weights], var = exposure, fam = svyglm_family)
impw_prs_adj_est <- .svyglm_pool_obese2(f_cov, imp_prs, weights = w[, weights], var = exposure, fam = svyglm_family)


ob2_t2d_plot <- tibble(
    mdm    = c(rep("Complete case", 4), rep ("Imputed", 4), rep("PRS-informed", 4)), # missing data method
    weight = c(rep(c("Unweighted", "Unweighted", "Weighted", "Weighted"), 3)),
    adjust = c(rep(c("Unadjusted", "Adjusted"), 6)),
    est    = exp(c(
        cc_un_est$est, cc_adj_est$est, ccw_un_est$est, ccw_adj_est$est,
        imp_un_est$est, imp_adj_est$est, impw_un_est$est, impw_adj_est$est,
        imp_prs_un_est$est, imp_prs_adj_est$est, impw_prs_un_est$est, impw_prs_adj_est$est
    )),
    lo    = exp(c(
        cc_un_est$lo, cc_adj_est$lo, ccw_un_est$lo, ccw_adj_est$lo,
        imp_un_est$lo, imp_adj_est$lo, impw_un_est$lo, impw_adj_est$lo,
        imp_prs_un_est$lo, imp_prs_adj_est$lo, impw_prs_un_est$lo, impw_prs_adj_est$lo
    )),
    hi    = exp(c(
        cc_un_est$hi, cc_adj_est$hi, ccw_un_est$hi, ccw_adj_est$hi,
        imp_un_est$hi, imp_adj_est$hi, impw_un_est$hi, impw_adj_est$hi,
        imp_prs_un_est$hi, imp_prs_adj_est$hi, impw_prs_un_est$hi, impw_prs_adj_est$hi
    ))
) |>
    ggplot(aes(x = mdm, y = est, color = weight, shape = weight)) +
    geom_rect(
        aes(xmin = -Inf, xmax = Inf, ymin = 3.4, ymax = 2.8),
        fill = "grey",
        alpha = 0.1,
        inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 1) +
    geom_pointrange(
        aes(ymin = lo, ymax = hi),
        position = position_dodge(width = 0.2),
        size = 1,
        linewidth = 1
    ) +
    # geom_point(position = position_dodge(width = 0.2), size = 3) +
    coord_flip() +
    facet_wrap(~ factor(adjust, c("Unadjusted", "Adjusted"))) +
    labs(
        title = "Type 2 diabetes coefficient for class-1 obesity by adjustment and missing data method",
        x = "Missing data method",
        y = "Odds ratio (95% CI)",
        caption = "Shaded region represents (covariate-adjusted) estimated range reported by Ganz et al. (2014)"
    ) +
    ms::scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.title = element_blank()
    )

# stack plots
patched <- (gluc_bmi_plot +
    labs(title = "A. BMI coefficient for glucose")) + 
    (t2d_bmi_plot +
        labs(title = "B. BMI coefficient for type 2 diabetes") +
        theme(legend.position = "none")) + 
    (ob1_t2d_plot +
        labs(title = "C. Type 2 diabetes coefficient for class-1 obesity ([30,35))") +
        theme(legend.position = "none")) +
    (ob2_t2d_plot +
        labs(title = "D. Type 2 diabetes coefficient for class-2 obesity ([35,40))") +
        theme(legend.position = "none")) +
        plot_layout(ncol = 1) +
        theme(legend.position = "none")

ggsave(
    plot = patched,
    filename = "results/stacked_real_world_comp.pdf",
    width = 10, height = 10,
    device = cairo_pdf
)
