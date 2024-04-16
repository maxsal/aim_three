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
make_missing <- function(data, intercept = -2.5, target_var, vars = NULL, target_var_wgt = NULL, var_wgts = NULL, scale_cont = TRUE, mech = "MAR") {
    dat      <- data.table::copy(data)
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
    dat
}

simul_diag <- function(cleaned, truth) {
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
