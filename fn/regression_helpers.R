# return nice results from glm estimates
.glm <- function(f, data, var, fam = "gaussian") {
    mod <- glm(f, data = data, family = fam)
    est <- mod$coefficients[[var]]
    ci <- suppressMessages({
        confint(mod)
    })[var, ]
    data.table(
        variable = var,
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
        variable = var,
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
        variable = var,
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
        variable = var,
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

    data.table(variable = var, est = pooled_mean, lo = lo, hi = hi, var = pooled_var)
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

    data.table(variable = var, est = pooled_mean, lo = lo, hi = hi, var = pooled_var)
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

    data.table(variance = var, est = pooled_mean, lo = lo, hi = hi, var = pooled_var)
}

# Custom function to apply .glm_pool with progress
map_glm_pool <- function(data_list, formula, exposure_variable) {
    map_dfr(
        data_list,
        \(x) .glm_pool(formula, x, exposure_variable),
        .progress = TRUE
    )
}

future_map_glm_pool <- function(data_list, formula, exposure_variable) {
    future_map_dfr(
        data_list,
        \(x) .glm_pool(formula, x, exposure_variable),
        .progress = TRUE,
        .options = furrr_options(seed = TRUE)
    )
}

map_svyglm_pool <- function(imp_list, formula, weights_list, var, fam = "gaussian") {
    map_dfr(
        seq_along(imp_list),
        \(x) .svyglm_pool(formula, imp_list[[x]], weights_list[[x]], var, fam),
        .progress = TRUE
    )
}

map_glm <- function(data_list, formula, exposure_variable) {
    map_dfr(
        data_list,
        \(x) .glm(formula, x, exposure_variable),
        .progress = TRUE
    )
}

future_map_glm <- function(data_list, formula, exposure_variable) {
    future_map_dfr(
        data_list,
        \(x) .glm(formula, x, exposure_variable),
        .progress = TRUE,
        .options = furrr_options(seed = TRUE)
    )
}

map_svyglm <- function(data_list, formula, exposure_variable) {
    map_dfr(
        data_list,
        \(x) .glm(formula, x, exposure_variable),
        .progress = TRUE
    )
}

monte_carlo_ci <- function(means, iterations, conf_level = 0.95) {
    # Sort the simulated means
    simulation_means <- sort(means)

    # Find the indices for the lower and upper bounds
    alpha <- (1 - conf_level) / 2
    lower_index <- floor(iterations * alpha)
    upper_index <- ceiling(iterations * (1 - alpha))

    # Extract the lower and upper bounds
    lower_bound <- simulation_means[max(1, lower_index)]
    upper_bound <- simulation_means[min(iterations, upper_index)]

    # Return the confidence interval
    data.table(
        est = mean(means),
        lo = lower_bound,
        hi = upper_bound
    )
}