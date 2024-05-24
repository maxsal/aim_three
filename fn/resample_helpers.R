# pull a resample of the data
re_sample <- function(x, size, prob = NULL, reps = 1, replace = TRUE, boot = FALSE) {
    if (boot) {
        size <- nrow(x)
        message(paste0("Bootstrapping with ", size, " samples (ignoring `size`` argument)"))
    }
    purrr::map(
        seq_len(reps),
        \(i) {
            index <- sample(x = nrow(x), size = size, replace = replace, prob = prob)
            x[index, ]
        },
        .progress = TRUE
    )
}



# make data missing at random
make_missing <- function(data, intercept = -2.5, target_var, vars = NULL, target_var_wgt = NULL, var_wgts = NULL, scale_cont = TRUE, mech = "MAR", mcar_prob = 0.25) {
    dat      <- data.table::copy(data)
    if (mech == "MCAR") {
        miss_ind <- rbinom(nrow(data), 1, prob = mcar_prob)
        dat[miss_ind == 1, target_var] <- NA
    } else {
        if (is.null(vars)) {
            vars <- setdiff(names(dat), target_var)
        }
        if (scale_cont) {
            num_vars  <- names(which(map_lgl(dat |> select(vars), is.numeric)))
            bin_vars  <- names(which(map_lgl(dat |> select(vars), \(x) all(na.omit(x) %in% 0:1))))
            cont_vars <- setdiff(num_vars, bin_vars)
            dat <- dat |>
                dplyr::mutate(dplyr::across(tidyselect::all_of(cont_vars), scale))
        }
        if (is.null(var_wgts)) {
            var_wgts <- rep(1, length(vars))
        }
        if (mech == "MAR") {
            score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(vars)) |> as.matrix())) %*% c(intercept, var_wgts)
        } else if (mech %in% c("MNAR", "NMAR")) {
            if (is.null(target_var_wgt)) {
                target_var_wgt <- 1
            }
            score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(c(vars, target_var))) |> as.matrix())) %*% c(intercept, var_wgts, target_var_wgt)
        }
        
        # score                   <- intercept + scale(x[["age"]])[, 1] + scale(x[["glucose"]])[, 1] + x[["female"]] + x[["smoke"]]
        miss_prob                      <- plogis(score)
        miss_ind                       <- rbinom(length(miss_prob), 1, miss_prob)
        dat[miss_ind == 1, target_var] <- NA
    }
    dat
}

make_ex_out_mar <- function(
    data,
    outcome,
    exposure,
    covs,
    ex_int = -3.28,
    out_int = -3.82) {
    dat <- data.table::copy(data)
    ex_vars <- c(outcome, covs)
    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(ex_vars)) |> as.matrix())) %*% c(ex_int, rep_len(1, length(ex_vars)))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_vars <- c(exposure, covs)
    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(out_vars)) |> as.matrix())) %*% c(out_int, rep_len(1, length(out_vars)))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

make_ex_out_mar2 <- function(
    data,
    outcome,
    exposure,
    covs,
    ex_int = -1.7,
    out_int = -1.9) {
    dat <- data.table::copy(data)
    ex_vars <- c(outcome, covs)
    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(ex_vars)) |> as.matrix())) %*% c(ex_int, c(1, 0.5, -0.25, 0.5, 0.25))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_vars <- c(exposure, covs)
    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(out_vars)) |> as.matrix())) %*% c(out_int, c(1, -0.5, 0.5, -0.25, -0.25))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

make_ex_out_mar3 <- function(data,
                                outcome,
                                exposure,
                                covs,
                                ex_int = -2.98,
                                out_int = -2.98) {
    dat <- data.table::copy(data)

    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(covs)) |> as.matrix())) %*% c(c(ex_int, rep_len(1, length(covs))))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(covs)) |> as.matrix())) %*% c(c(out_int, rep_len(1, length(covs))))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

make_ex_out_mar4 <- function(data,
                             outcome,
                             exposure,
                             covs,
                             ex_int = -2.98,
                             out_int = -2.98) {
    dat <- data.table::copy(data)

    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(covs)) |> as.matrix())) %*% c(c(ex_int, rep_len(1, length(covs))))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_vars <- c(exposure, covs)
    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(out_vars)) |> as.matrix())) %*% c(c(out_int, rep_len(1, length(out_vars))))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

make_ex_out_mar5 <- function(data,
         outcome,
         exposure,
         covs,
         ex_int = -3.21,
         out_int = -2.98
) {
    dat <- data.table::copy(data)

    ex_vars <- c(outcome, covs)
    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(ex_vars)) |> as.matrix())) %*% c(c(ex_int, rep_len(1, length(ex_vars))))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(covs)) |> as.matrix())) %*% c(c(out_int, rep_len(1, length(covs))))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

make_ex_out_mnar <- function(
    data,
    outcome,
    exposure,
    covs,
    ex_int = -1.78,
    out_int = -2.68) {
    dat <- data.table::copy(data)

    ex_vars <- c(exposure, outcome, covs)
    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(ex_vars)) |> as.matrix())) %*% c(ex_int, 1, 1, 1, -0.25, -0.25, -0.25)
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_vars <- c(outcome, exposure, covs)
    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(out_vars)) |> as.matrix())) %*% c(out_int, 1, 1, 0.5, 0.25, 0.5, 0.25)
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

make_ex_out_mnar3 <- function(
    data,
    outcome,
    exposure,
    covs,
    ex_int = -3.96,
    out_int = -3.96) {
    dat <- data.table::copy(data)

    ex_vars <- c(outcome, exposure, covs)
    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(ex_vars)) |> as.matrix())) %*% c(ex_int, rep_len(1, length(ex_vars)))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_vars <- c(exposure, outcome, covs)
    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(out_vars)) |> as.matrix())) %*% c(out_int, rep_len(1, length(out_vars)))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}


make_ex_out_mnar5 <- function(
    data,
    outcome,
    exposure,
    covs,
    ex_int = -3.95,
    out_int = -3.22) {
    dat <- data.table::copy(data)

    ex_vars <- c(outcome, exposure, covs)
    ex_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(ex_vars)) |> as.matrix())) %*% c(ex_int, rep_len(1, length(ex_vars)))
    ex_miss_prob <- plogis(ex_score)
    ex_miss_ind <- rbinom(length(ex_miss_prob), 1, ex_miss_prob)

    out_vars <- c(outcome, covs)
    out_score <- as.matrix(cbind(1, dat |> dplyr::select(tidyselect::all_of(out_vars)) |> as.matrix())) %*% c(out_int, rep_len(1, length(out_vars)))
    out_miss_prob <- plogis(out_score)
    out_miss_ind <- rbinom(length(out_miss_prob), 1, out_miss_prob)

    dat[ex_miss_ind == 1, exposure] <- NA
    dat[out_miss_ind == 1, outcome] <- NA

    return(dat)
}

# simul_diag <- function(cleaned, truth) {
#     raw_bias <- mean(cleaned$est - truth)
#     per_bias <- 100 * abs(raw_bias / truth)
#     cov_rate <- mean(cleaned$lo < truth & cleaned$hi > truth)
#     avg_wide <- mean(cleaned$hi - cleaned$lo)
#     rmse <- sqrt(mean((cleaned$est - truth)^2))
#     data.table(
#         raw_bias = raw_bias,
#         per_bias = per_bias,
#         cov_rate = cov_rate,
#         avg_wide = avg_wide,
#         rmse     = rmse
#     )
# }

  simul_diag <- function(cleaned, truth) {
      raw_bias <- mean(cleaned$est - truth)
      per_bias <- mean(100 * abs(raw_bias / truth))
      cov_rate <- mean(cleaned$lo < truth & cleaned$hi > truth)
      avg_wide <- mean(cleaned$hi - cleaned$lo)
      rmse <- sqrt(mean((cleaned$est - truth)^2))
      data.table(
          adj       = cleaned[, unique(adj)],
          miss_data = cleaned[, unique(miss_data)],
          weight    = cleaned[, unique(weight)],
          size      = cleaned[, unique(size)],
          method    = cleaned[, unique(method)],
          raw_bias  = raw_bias,
          per_bias  = per_bias,
          cov_rate  = cov_rate,
          avg_wide  = avg_wide,
          rmse      = rmse
      )
  }

    simul_diag2 <- function(cleaned, truth) {
        raw_bias <- mean(cleaned$est - truth)
        per_bias <- mean(100 * abs(raw_bias / truth))
        cov_rate <- mean(cleaned$lo < truth & cleaned$hi > truth)
        avg_wide <- mean(cleaned$hi - cleaned$lo)
        rmse <- sqrt(mean((cleaned$est - truth)^2))
        data.table(
            raw_bias  = raw_bias,
            per_bias  = per_bias,
            cov_rate  = cov_rate,
            avg_wide  = avg_wide,
            rmse      = rmse
        )
    }