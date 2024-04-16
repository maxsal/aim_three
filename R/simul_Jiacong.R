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

data <- d |> select(age, female, glucose, bmi, smoke)

scaled <- data |>
    mutate(
        age = scale(age)[, 1],
        glucose = scale(glucose)[, 1],
        bmi = scale(bmi)[, 1]
    )

fit <- lm(data = scaled, formula = glucose ~ bmi)
sig <- sqrt(var(fit$residuals))
hist(fit$residuals)

hist(scaled$glucose)
hist(predict(fit, data = scaled) + rnorm(nrow(scaled), 0, sig))

quick_eval <- function(cleaned, truth) {
    raw_bias <- mean(cleaned$est - truth)
    per_bias <- 100 * abs(raw_bias / truth)
    cov_rate <- mean(cleaned$lo < truth & cleaned$hi > truth)
    avg_wide <- mean(cleaned$hi - cleaned$lo)
    rmse <- sqrt(mean((cleaned$est - truth)^2))
    data.table(
        raw_bias = raw_bias,
        per_bias = per_bias,
        cov_rate = cov_rate,
        avg_wide = avg_wide,
        rmse     = rmse
    )
}

######################### Approach 2 #########################
rho         <- lm(data = scaled, formula = glucose ~ bmi)$coefficients[["bmi"]]
reps        <- 100
cover       <- c()
cc_cover    <- c()
data_list   <- list()
imps        <- list()
miss_prop   <- c()
cc_results  <- list()
imp_results <- list()
exposure    <- "bmi"

# sample <- scaled
for (rep in 1:reps) {
    # step 1: set.seed
    set.seed(rep)
    # step 2: generate outcome
    glucose_new <- predict(fit, data = scaled) + rnorm(nrow(scaled), 0, sig)
    # # step 3: affix outcome to data
    sample <- scaled |>
        mutate(glucose = glucose_new)
    # step 4: generate missing data
    scores                        <- -1 - sample[["glucose"]] - sample[["age"]] - sample[["female"]]
    miss_probs                    <- plogis(scores)
    miss_ind                      <- rbinom(length(miss_probs), 1, miss_probs)
    man                           <- sample
    man[miss_ind == 1, "bmi"] <- NA
    # step 5: impute missing data
    imp <- mice(man, m = 5, method = "pmm", print = FALSE)
    # step 6: get cc and imp estimates
        ## cc
            cc_fit      <- lm(glucose ~ bmi, data = man)
            cc_result_i <- summary(cc_fit)$coefficients[2, c("Estimate", "Std. Error")]
            cc_ci_low   <- cc_result_i[1] - qnorm(0.975) * cc_result_i[2]
            cc_ci_up    <- cc_result_i[1] + qnorm(0.975) * cc_result_i[2]
            cc_cover_i  <- (cc_ci_low < rho) & (cc_ci_up > rho)
        ## imp
            pool_fit <- pool(with(imp, lm(glucose ~ bmi)))
            result_i <- summary(pool_fit)[2, c("estimate", "std.error")]
            ci_low   <- result_i[[1]] - 2 * result_i[[2]]
            ci_up    <- result_i[[1]] + 2 * result_i[[2]]
            cover_i  <- (ci_low < rho) & (ci_up > rho)
    # step 7: save results
        # - data
        # - imp object
        # - average missing proportion
        # - beta estimates with confidence intervals
        # save the results
        data_list[[rep]]   <- man
        imps[[rep]]        <- imp
        miss_prop         <- c(miss_prop, mean(miss_ind))
        cc_results[[rep]]  <- data.table(est = cc_result_i[1], lo = cc_ci_low, hi = cc_ci_up)
        imp_results[[rep]] <- data.table(est = result_i[[1]], lo = ci_low, hi = ci_up)

        cc_cover <- c(cc_cover, cc_cover_i)
        cover <- c(cover, cover_i)
        print(paste0(
            "I: ", rep,
            "; CCcoverage: ", trimws(format(round(mean(cc_cover), 3), nsmall = 3)),
            "; IMPcoverage: ", trimws(format(round(mean(cover), 3), nsmall = 3))
            ))

}

cleaned_cc  <- rbindlist(cc_results)
cleaned_imp <- rbindlist(imp_results)

bind_rows(
    cbind(data.table(analysis = "Complete case"), quick_eval(cleaned_cc, rho)),
    cbind(data.table(analysis = "Imputed"), quick_eval(cleaned_imp, rho))
)

tibble(
    index = 1:reps,
    coverage = dplyr::cummean(cc_cover)
) |>
    ggplot(aes(x = index, y = coverage)) +
    geom_line() +
    ms::theme_ms()

bind_rows(
    cbind(data.table(analysis = "Complete case"), cleaned_cc),
    cbind(data.table(analysis = "Imputed"), cleaned_imp)
) |>
    ggplot(aes(x = est, color = analysis)) +
    geom_vline(xintercept = rho, linetype = "dashed", linewidth = 1) +
    geom_density(linewidth = 1) +
    ms::theme_ms()
ggsave(
    filename = "bin/test_plot.pdf",
    width = 6, height = 4,
    device = cairo_pdf
)
#### TRY WITH AMPUTE ####
rho <- lm(data = scaled, formula = glucose ~ bmi)$coefficients[["bmi"]]
reps <- 100
cover <- c()
cc_cover <- c()
data_list <- list()
imps <- list()
miss_prop <- c()
cc_results <- list()
imp_results <- list()
exposure <- "bmi"

# tmp_scaled <- scaled
# tmp_index <- sample(1:nrow(tmp_scaled), 1000, replace = FALSE)
# scaled <- scaled[tmp_index, ]
# test <- ampute(scaled, prop = 0.5, mech = "MAAR", patterns = c(1, 1, 0, 1, 1))

for (rep in 1:reps) {
    # step 1: set.seed
    set.seed(rep)
    # step 2: ampute data
    man <- ampute(scaled, prop = 0.5, mech = "MNAR", freq = c(0.1, 0.1, 0.3, 0.4, 0.1))$amp
    # step 5: impute missing data
    imp <- mice(man, m = 5, method = "pmm", print = FALSE)
    # step 6: get cc and imp estimates
    ## cc
    cc_fit <- lm(glucose ~ bmi, data = man)
    cc_result_i <- summary(cc_fit)$coefficients[2, c("Estimate", "Std. Error")]
    cc_ci_low <- cc_result_i[1] - qnorm(0.975) * cc_result_i[2]
    cc_ci_up <- cc_result_i[1] + qnorm(0.975) * cc_result_i[2]
    cc_cover_i <- (cc_ci_low < rho) & (cc_ci_up > rho)
    ## imp
    pool_fit <- pool(with(imp, lm(glucose ~ bmi)))
    result_i <- summary(pool_fit)[2, c("estimate", "std.error")]
    ci_low <- result_i[[1]] - 2 * result_i[[2]]
    ci_up <- result_i[[1]] + 2 * result_i[[2]]
    cover_i <- (ci_low < rho) & (ci_up > rho)
    # step 7: save results
    # - data
    # - imp object
    # - average missing proportion
    # - beta estimates with confidence intervals
    # save the results
    data_list[[rep]] <- man
    imps[[rep]] <- imp
    miss_prop <- c(miss_prop, mean(miss_ind))
    cc_results[[rep]] <- data.table(est = cc_result_i[1], lo = cc_ci_low, hi = cc_ci_up)
    imp_results[[rep]] <- data.table(est = result_i[[1]], lo = ci_low, hi = ci_up)
    results <- rbind(results, result_i)

    cc_cover <- c(cc_cover, cc_cover_i)
    cover <- c(cover, cover_i)
    print(paste0(
        "I:", rep,
        ": CCcoverage: ", trimws(format(round(mean(cc_cover), 3), nsmall = 3)),
        "; IMPcoverage: ", trimws(format(round(mean(cover), 3), nsmall = 3))
    ))
}

cleaned_cc <- rbindlist(cc_results)
cleaned_imp <- rbindlist(imp_results)

bind_rows(
    cbind(data.table(analysis = "Complete case"), quick_eval(cleaned_cc, rho)),
    cbind(data.table(analysis = "Imputed"), quick_eval(cleaned_imp, rho))
)

# ####
# for (rep in 1:reps) {
#     # generate the data

#     set.seed(rep)
#     # sample <- MASS::mvrnorm(
#     #     n = 1000,
#     #     mu = c(0, 0, 0),
#     #     Sigma = matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), ncol = 3)
#     # ) |> as.data.frame()
#     # names(sample) <- c("Y1", "Y2", "X1")
#     glucose_new <- predict(fit, data = scaled) + rnorm(nrow(scaled), 0, sig)
#     sample <- scaled |>
#         mutate(glucose = glucose_new)

#     # summary(lm(Y1 ~ X1, data = sample)) # adjusted R2 is <0.05, a bit small even though X1 is a strong predictor of Y1
#     # summary(lm(Y1 ~ X1+Y2, data = sample))

#     # generate missing data
#     scores     <- -1 - sample[["bmi"]] - sample[["age"]] - sample[["female"]]
#     miss_probs <- plogis(scores)
#     miss_ind   <- rbinom(length(miss_probs), 1, miss_probs)

#     man <- sample
#     man[miss_ind == 1, "glucose"] <- NA

#     # imputed data using mice and trace=0
#     imp <- mice(man, m = 5, method = "pmm", print = FALSE)
#     pool_fit <- pool(with(imp, lm(glucose ~ bmi)))
#     result_i <- summary(pool_fit)[2, c("estimate", "std.error")]
#     ci_low <- result_i[1] - 2 * result_i[2]
#     ci_up <- result_i[1] + 2 * result_i[2]
#     cover_i <- (ci_low < rho) & (ci_up > rho)

#     # save the results
#     results <- rbind(results, result_i)

#     cover <- c(cover, cover_i)
#     print(paste0(rep, "-th run: ", mean(cover)))
# }
