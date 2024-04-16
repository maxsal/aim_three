ms::libri(ms, qs, data.table, tidyverse, glue, survey, mice)
source("fn/regression_helpers.R")

# data -------------------------------------------------------------------------
## read data
d <- qread("data/private/mgi_demo_comb_20230322.qs")
w <- qread("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/mgi/20220822/weightsx_20220822_comb.qs")
p <- qread("data/private/mgi_pim0x_20230322.qs")
l <- qread("data/private/mgi_labs_cleaned_20230322.qs")

prs <- qread("data/private/ex_prs_processed.qs")
prs <- prs[, .(id, bmi_prs = bmi, t2d_prs = t2d, glucose_prs = glucose)]

## remove people with missing weight values
w <- w[!is.na(ip_selection), ]

## identify the intersection
ids <- reduce(list(d[, id], w[, id], p[, id]), intersect)

d <- d[id %in% ids, ]
w <- w[id %in% ids, ]

## merge data
comb <- reduce(
    list(
        d[, .(id, age, female, nhw = as.numeric(race_eth == "NHW"), smoke = smoke_ever, bmi = bmi_med)],
        w[, .(id, weight = ip_selection)],
        p[, .(id, t2d = EM_202.2)],
        l[, .(id, glucose = glucose_med)]
    ),
    left_join,
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

glu_bmi_p_nhanes_est <- tibble(
    adj = factor(c("Unadjusted", "Adjusted"), c("Unadjusted", "Adjusted")),
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

res_tab |>
    mutate(
        weight = if_else(weight == "Yes", "Weighted", "Unweighted")
    ) |>
    # filter(mdm != "Imputation w/ PRS") |>
    ggplot(aes(x = mdm, y = est, ymin = lo, ymax = hi, color = weight)) +
        geom_rect(
        data = glu_bmi_p_nhanes_est,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = NULL, y = NULL),
        fill = "grey",
        alpha = 0.5,
        show.legend = FALSE,
        inherit.aes = FALSE
    ) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
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
    "bin/bmi_glucose_coef_w_range.pdf",
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

tmp <- glm(f_cov, data = dat)

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
    ggplot(aes(x = mdm, y = est, ymin = lo, ymax = hi, color = weight)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
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
    "bin/glucose_bmi_coef.pdf",
    width = 6, height = 5,
    device = cairo_pdf
)
