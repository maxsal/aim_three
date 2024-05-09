ms::libri(ms, qs, data.table, tidyverse, glue, survey, mice)
source("fn/regression_helpers.R")

# data -------------------------------------------------------------------------
## read data
d <- qread("data/private/mgi_demo_comb_20230322.qs")
w <- qread("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/mgi/20220822/weightsx_20220822_comb.qs")
p <- qread("data/private/mgi_pim0x_20230322.qs")
l <- qread("data/private/mgi_labs_cleaned_20230322.qs")

prs <- qread("data/private/ex_prs_processed.qs")
prs <- prs[, .(id, bmi_prs = bmi, glucose_prs = glucose)]

## remove people with missing weight values
w <- w[!is.na(ip_selection), ]

## keep people with non-missing weights
d   <- d[id %in% w[, id], ]
w   <- w[id %in% w[, id], ]
p   <- p[id %in% w[, id], ]
l   <- l[id %in% w[, id], ]
prs <- prs[id %in% w[, id], ]

## merge data
comb <- reduce(
    list(
        d[, .(id, age, female, nhw = as.numeric(race_eth == "NHW"), smoke = smoke_ever, bmi = bmi_med)],
        w[, .(id, weight = ip_selection)],
        p[, .(id, t2d = EM_202.2)],
        l[, .(id, glucose = glucose_med)]
    ),
    full_join,
    by = "id"
)

lapply(comb, \(x) sum(is.na(x)))

# Glucose ~ BMI analysis -------------------------------------------------------
outcome    <- "glucose"
exposure   <- "bmi"
covariates <- c("age", "female", "nhw", "smoke")
dat <- comb |>
    select(all_of(c(outcome, exposure, covariates)))
non_outcome <- names(dat) |> setdiff(outcome)
f     <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(non_outcome, collapse = " + "))

# COMPLETE CASE
# standard
cc_un_est  <- .glm(f, data = dat, var = exposure)
cc_adj_est <- .glm(f_cov, data = dat, var = exposure)

# with weights
datw <- cbind(comb[, .(id)], dat) |>
    left_join(comb[, .(id, weight)], by = "id")
cc_dsn <- svydesign(
    ids     = ~1,
    data    = datw,
    weights = ~weight
)

ccw_un_est  <- .svyglm(f, var = exposure, design = cc_dsn)
ccw_adj_est <- .svyglm(f_cov, var = exposure, design = cc_dsn)

# tmp <- glm(f_cov, data = dat)

# IMPUTATION
imp <- mice::mice(dat, m = 5, seed = 123, printFlag = FALSE)

imp_un_est  <- .glm_pool(f = f, imp = imp, var = exposure)
imp_adj_est <- .glm_pool(f = f_cov, imp = imp, var = exposure)

impw_un_est  <- .svyglm_pool(f, imp, weights = w[, ip_selection], var = exposure)
impw_adj_est <- .svyglm_pool(f_cov, imp, weights = w[, ip_selection], var = exposure)

# IMPUTATION W PRS
dat_prs <- left_join(
    cbind(comb[, .(id)], dat),
    prs,
    by = "id"
)[, !c("id", "t2d_prs")] 

impp <- mice::mice(dat_prs, m = 5, seed = 123, printFlag = FALSE)

impp_un_est  <- .glm_pool(f = f, imp = impp, var = exposure)
impp_adj_est <- .glm_pool(f = f_cov, imp = impp, var = exposure)

imppw_un_est  <- .svyglm_pool(f, impp, weights = w[, ip_selection], var = exposure)
imppw_adj_est <- .svyglm_pool(f_cov, impp, weights = w[, ip_selection], var = exposure)

# remove people with missing PRS
dat_prs_sub    <- cbind(comb[, .(id)], dat_prs)[!(is.na(bmi_prs) | is.na(glucose_prs)), ]
dat_prs_sub    <- left_join(dat_prs_sub, w[, .(id, weights = ip_selection)], by = "id")
dat_prs_sub_id <- dat_prs_sub$id
dat_prs_sub_w  <- dat_prs_sub$weights
dat_prs_sub    <- dat_prs_sub[, !c("id", "weights")]

imps <- mice::mice(dat_prs_sub, m = 5, seed = 123, printFlag = FALSE)

imps_un_est  <- .glm_pool(f = f, imp = imps, var = exposure)
imps_adj_est <- .glm_pool(f = f_cov, imp = imps, var = exposure)

impsw_un_est  <- .svyglm_pool(f, imps, weights = dat_prs_sub_w, var = exposure)
impsw_adj_est <- .svyglm_pool(f_cov, imps, weights = dat_prs_sub_w, var = exposure)

res_tab <- tribble(
    ~mdm, ~adj, ~weight, ~est, ~lo, ~hi,
    "Complete case", "Unadjusted", "No", cc_un_est$est, cc_un_est$lo, cc_un_est$hi,
    "Complete case", "Adjusted", "No", cc_adj_est$est, cc_adj_est$lo, cc_adj_est$hi,
    "Complete case", "Unadjusted", "Yes", ccw_un_est$est, ccw_un_est$lo, ccw_un_est$hi,
    "Complete case", "Adjusted", "Yes", ccw_adj_est$est, ccw_adj_est$lo, ccw_adj_est$hi,
    "Imputation", "Unadjusted", "No", imp_un_est$est, imp_un_est$lo, imp_un_est$hi,
    "Imputation", "Adjusted", "No", imp_adj_est$est, imp_adj_est$lo, imp_adj_est$hi,
    "Imputation", "Unadjusted", "Yes", impw_un_est$est, impw_un_est$lo, impw_un_est$hi,
    "Imputation", "Adjusted", "Yes", impw_adj_est$est, impw_adj_est$lo, impw_adj_est$hi,
    "Imputation w/ PRS", "Unadjusted", "No", impp_un_est$est, impp_un_est$lo, impp_un_est$hi,
    "Imputation w/ PRS", "Adjusted", "No", impp_adj_est$est, impp_adj_est$lo, impp_adj_est$hi,
    "Imputation w/ PRS", "Unadjusted", "Yes", imppw_un_est$est, imppw_un_est$lo, imppw_un_est$hi,
    "Imputation w/ PRS", "Adjusted", "Yes", imppw_adj_est$est, imppw_adj_est$lo, imppw_adj_est$hi,
    "Imputation w/ PRS (subset)", "Unadjusted", "No", imps_un_est$est, imps_un_est$lo, imps_un_est$hi,
    "Imputation w/ PRS (subset)", "Adjusted", "No", imps_adj_est$est, imps_adj_est$lo, imps_adj_est$hi,
    "Imputation w/ PRS (subset)", "Unadjusted", "Yes", impsw_un_est$est, impsw_un_est$lo, impsw_un_est$hi,
    "Imputation w/ PRS (subset)", "Adjusted", "Yes", impsw_adj_est$est, impsw_adj_est$lo, impsw_adj_est$hi
) |>
    mutate(
        adj = factor(adj, levels = c("Unadjusted", "Adjusted")),
        weight = factor(weight, levels = c("No", "Yes")),
        mdm = factor(mdm, levels = c("Complete case", "Imputation", "Imputation w/ PRS", "Imputation w/ PRS (subset)"))
    )

# glu_bmi_p_nhanes_est <- tibble(
#     adj = factor(c("Unadjusted", "Adjusted"), c("Unadjusted", "Adjusted")),
#     xmin = c(-Inf, -Inf),
#     xmax = c(Inf, Inf),
#     ymin = c(0.624, 0.601),
#     ymax = c(1.098, 1.061),
#     est = c(NA_real_, NA_real_),
#     lo = c(NA_real_, NA_real_),
#     hi = c(NA_real_, NA_real_),
#     mdm = c(NA_character_, NA_character_),
#     weight = c(NA_character_, NA_character_)
# )

res_tab |>
    mutate(
        weight = if_else(weight == "Yes", "Weighted", "Unweighted"),
        mdm = case_when(
            mdm == "Complete case" ~ "Complete case",
            mdm == "Imputation" ~ "OD-imputed",
            mdm == "Imputation w/ PRS" ~ "OD-PRS-imputed",
            mdm == "Imputation w/ PRS (subset)" ~ "OD-PRS-imputed (subset)"
        ),
        adj = factor(case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        ), c("Unadjusted", "Covariate-adjusted"))
    ) |>
    # filter(mdm != "Imputation w/ PRS") |>
    ggplot(aes(x = mdm, y = est, ymin = lo, ymax = hi, color = weight, shape = weight)) +
    # geom_rect(
    #     data = glu_bmi_p_nhanes_est,
    #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = NULL, y = NULL),
    #     fill = "grey",
    #     alpha = 0.5,
    #     show.legend = FALSE,
    #     inherit.aes = FALSE
    # ) +
    geom_pointrange(size = 1, linewidth = 1, position = position_dodge(width = 0.5)) +
        labs(
            title = "BMI coefficient for glucose",
            x = "Missing data method",
            y = "Beta coefficient (95% CI)"
            # caption = "Shaded region represents estimated range using NHANES 2017-2020 (P) data (fasting glucose)"
        ) +
    # ylim(0.4, 0.8) +
    facet_grid(~adj) +
    scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(
    "bin/glucose_bmi_coef.pdf",
    width = 6, height = 5,
    device = cairo_pdf
)

# BMI ~ glucose analysis -------------------------------------------------------
outcome <- "bmi"
exposure <- "glucose"
covariates <- c("age", "female", "nhw", "smoke")
dat <- comb |>
    select(all_of(c(outcome, exposure, covariates)))
non_outcome <- names(dat) |> setdiff(outcome)
f <- paste0(outcome, " ~ ", exposure)
f_cov <- paste0(outcome, " ~ ", paste(non_outcome, collapse = " + "))

# COMPLETE CASE
# standard
cc_un_est <- .glm(f, data = dat, var = exposure)
cc_adj_est <- .glm(f_cov, data = dat, var = exposure)

# with weights
datw <- cbind(comb[, .(id)], dat) |>
    left_join(comb[, .(id, weight)], by = "id")
cc_dsn <- svydesign(
    ids     = ~1,
    data    = datw,
    weights = ~weight
)

ccw_un_est <- .svyglm(f, var = exposure, design = cc_dsn)
ccw_adj_est <- .svyglm(f_cov, var = exposure, design = cc_dsn)

# tmp <- glm(f_cov, data = dat)

# IMPUTATION
imp <- mice::mice(dat, m = 5, seed = 123, printFlag = FALSE)

imp_un_est <- .glm_pool(f = f, imp = imp, var = exposure)
imp_adj_est <- .glm_pool(f = f_cov, imp = imp, var = exposure)

impw_un_est <- .svyglm_pool(f, imp, weights = w[, ip_selection], var = exposure)
impw_adj_est <- .svyglm_pool(f_cov, imp, weights = w[, ip_selection], var = exposure)

# IMPUTATION W PRS
dat_prs <- left_join(
    cbind(comb[, .(id)], dat),
    prs,
    by = "id"
)[, !c("id", "t2d_prs")]

impp <- mice::mice(dat_prs, m = 5, seed = 123, printFlag = FALSE)

impp_un_est <- .glm_pool(f = f, imp = impp, var = exposure)
impp_adj_est <- .glm_pool(f = f_cov, imp = impp, var = exposure)

imppw_un_est <- .svyglm_pool(f, impp, weights = w[, ip_selection], var = exposure)
imppw_adj_est <- .svyglm_pool(f_cov, impp, weights = w[, ip_selection], var = exposure)

# remove people with missing PRS
dat_prs_sub <- cbind(comb[, .(id)], dat_prs)[!(is.na(bmi_prs) | is.na(glucose_prs)), ]
dat_prs_sub <- left_join(dat_prs_sub, w[, .(id, weights = ip_selection)], by = "id")
dat_prs_sub_id <- dat_prs_sub$id
dat_prs_sub_w <- dat_prs_sub$weights
dat_prs_sub <- dat_prs_sub[, !c("id", "weights")]

imps <- mice::mice(dat_prs_sub, m = 5, seed = 123, printFlag = FALSE)

imps_un_est <- .glm_pool(f = f, imp = imps, var = exposure)
imps_adj_est <- .glm_pool(f = f_cov, imp = imps, var = exposure)

impsw_un_est <- .svyglm_pool(f, imps, weights = dat_prs_sub_w, var = exposure)
impsw_adj_est <- .svyglm_pool(f_cov, imps, weights = dat_prs_sub_w, var = exposure)

res_tab <- tribble(
    ~mdm, ~adj, ~weight, ~est, ~lo, ~hi,
    "Complete case", "Unadjusted", "No", cc_un_est$est, cc_un_est$lo, cc_un_est$hi,
    "Complete case", "Adjusted", "No", cc_adj_est$est, cc_adj_est$lo, cc_adj_est$hi,
    "Complete case", "Unadjusted", "Yes", ccw_un_est$est, ccw_un_est$lo, ccw_un_est$hi,
    "Complete case", "Adjusted", "Yes", ccw_adj_est$est, ccw_adj_est$lo, ccw_adj_est$hi,
    "Imputation", "Unadjusted", "No", imp_un_est$est, imp_un_est$lo, imp_un_est$hi,
    "Imputation", "Adjusted", "No", imp_adj_est$est, imp_adj_est$lo, imp_adj_est$hi,
    "Imputation", "Unadjusted", "Yes", impw_un_est$est, impw_un_est$lo, impw_un_est$hi,
    "Imputation", "Adjusted", "Yes", impw_adj_est$est, impw_adj_est$lo, impw_adj_est$hi,
    "Imputation w/ PRS", "Unadjusted", "No", impp_un_est$est, impp_un_est$lo, impp_un_est$hi,
    "Imputation w/ PRS", "Adjusted", "No", impp_adj_est$est, impp_adj_est$lo, impp_adj_est$hi,
    "Imputation w/ PRS", "Unadjusted", "Yes", imppw_un_est$est, imppw_un_est$lo, imppw_un_est$hi,
    "Imputation w/ PRS", "Adjusted", "Yes", imppw_adj_est$est, imppw_adj_est$lo, imppw_adj_est$hi,
    "Imputation w/ PRS (subset)", "Unadjusted", "No", imps_un_est$est, imps_un_est$lo, imps_un_est$hi,
    "Imputation w/ PRS (subset)", "Adjusted", "No", imps_adj_est$est, imps_adj_est$lo, imps_adj_est$hi,
    "Imputation w/ PRS (subset)", "Unadjusted", "Yes", impsw_un_est$est, impsw_un_est$lo, impsw_un_est$hi,
    "Imputation w/ PRS (subset)", "Adjusted", "Yes", impsw_adj_est$est, impsw_adj_est$lo, impsw_adj_est$hi
) |>
    mutate(
        adj = factor(adj, levels = c("Unadjusted", "Adjusted")),
        weight = factor(weight, levels = c("No", "Yes")),
        mdm = factor(mdm, levels = c("Complete case", "Imputation", "Imputation w/ PRS", "Imputation w/ PRS (subset)"))
    )

res_tab |>
    mutate(
        weight = if_else(weight == "Yes", "Weighted", "Unweighted")
    ) |>
    # filter(mdm != "Imputation w/ PRS") |>
    ggplot(aes(x = mdm, y = est, ymin = lo, ymax = hi, color = weight, shape = weight)) +
    geom_pointrange(size = 1, linewidth = 1, position = position_dodge(width = 0.5)) +
    labs(
        title = "Glucose coefficient for BMI",
        x = "Missing data method",
        y = "Beta coefficient (95% CI)"
    ) +
    facet_grid(~adj) +
    scale_color_ms() +
    ms::theme_ms() +
    theme(
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(
    "bin/bmi_glucose_coef.pdf",
    width = 6, height = 5,
    device = cairo_pdf
)

# MISSINGNESS PLOT AND TABLE ---------------------------------------------------
full <- dat_prs
full$missing <- ifelse(complete.cases(full), "Fully observed", "Any missing")
full$cc <- ifelse(complete.cases(full), "Complete case (adjusted)",
                  ifelse(complete.cases(full[, .(glucose, bmi)]), "Complete case (unadjusted)", NA))


# OVERALL
# Complete case (adjusted)
    # Any missing
    # No missing
# Complete case (adjusted) - with PRS

sub <- na.omit(full)

libri(gtsummary)

full |>
    mutate(
        `BMI PRS` = scale(bmi_prs)[, 1],
        `Glucose PRS` = scale(glucose_prs)[, 1]
    ) |>
    select(
        Age = age,
        Female = female,
        `Non-Hispanic White` = nhw,
        `Smoking status (ever)` = smoke,
        `BMI` = bmi,
        `Glucose` = glucose,
        `BMI PRS`,
        `Glucose PRS`, missing
    ) |>
    tbl_summary(
        by = missing,
        missing_text = "Missing",
        statistic = list(
            all_continuous() ~ "{mean} ({sd})",
            all_dichotomous() ~ "{p} ({n})"),
        digits = list(
            all_continuous() ~ c(1, 1),
            all_dichotomous() ~ c(1, 0),
            contains("PRS") ~ c(3, 3)
        )) |>
    add_p() |>
    add_overall()

full |>
    mutate(
        `BMI PRS` = scale(bmi_prs)[, 1],
        `Glucose PRS` = scale(glucose_prs)[, 1]
    ) |>
    drop_na(bmi_prs, glucose_prs) |>
    select(
        Age = age,
        Female = female,
        `Non-Hispanic White` = nhw,
        `Smoking status (ever)` = smoke,
        `BMI` = bmi,
        `Glucose` = glucose,
        `BMI PRS`,
        `Glucose PRS`, missing
    ) |>
    tbl_summary(
        by = missing,
        missing_text = "Missing",
        statistic = list(
            all_continuous() ~ "{mean} ({sd})",
            all_dichotomous() ~ "{p} ({n})"
        ),
        digits = list(
            all_continuous() ~ c(1, 1),
            all_dichotomous() ~ c(1, 0),
            contains("PRS") ~ c(3, 3)
        )
    ) |>
    add_p() |>
    add_overall()

full |>
    select(-missing) |>
    rename(
        `Glucose (lab)` = glucose,
        `BMI` = bmi,
        Age = age,
        Female = female,
        `Non-Hispanic White` = nhw,
        `Smoking status (ever)` = smoke,
        `Glucose PRS` = glucose_prs,
        `BMI PRS` = bmi_prs
    ) |>
    map_dbl(\(x) sum(is.na(x)) / length(x)) |>
    enframe() |>
    arrange(desc(value)) |>
    ggplot(aes(x = reorder(name, value), y = value)) +
    geom_segment(aes(xend = name, yend = 0), linewidth = 2, color = "grey") +
    geom_point(size = 4) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 0.5)) +
    labs(
        title = "Missingness of relavent variables in MGI",
        x = "Variable",
        y = "Missingness (%)"
    ) +
    coord_flip() +
    theme_ms() +
    theme(
        panel.grid.major.y = element_blank()
    )
ggsave(
    "results/missingness_plot.pdf",
    width = 6,
    height = 5,
    device = cairo_pdf
)
