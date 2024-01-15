# Functions for helping curate data

## Filtering out extreme values ------------------------------------------------
## Default
# approach from Ma,...,Fritsche 2022 (doi: 10.1016/j.ajhg.2022.09.001)
exprs_clean <- function(x) {
    x_qt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
    x_iqr <- IQR(x, na.rm = TRUE)
    x[(x > (x_qt[2] + 1.5*x_iqr))] <- NA
    x[(x < (x_qt[1] - 1.5*x_iqr))] <- NA
    x
}

## Alternative
# approach from Li, Chen & Moore 2019 (PMID: 31329892)
sd4_clean <- function(x) {
    mn <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    x[x > mn + 4 * sd] <- NA
    x[x < mn - 4 * sd] <- NA
    x
}

# helper function to read number of rows in a file
get_line_count <- function(file) {
    f <- file(file, open = "rb")
    nlines <- 0L
    while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
        nlines <- nlines + sum(chunk == as.raw(10L))
    }
    nlines
}


## Cleaning the lab data based on LOINC code -----------------------------------
clean_labs <- function(lab_data, loinc_code, name) {
    data_for_cleaning <- lab_data[loinc == loinc_code, ]

    # starting counts
    suppressWarnings({
        n_total_obs <- data_for_cleaning[, .N]
        n_total_ppl <- data_for_cleaning[, .N, id][, .N]
        n_miss_num_conv_obs <- sum(is.na(as.numeric(data_for_cleaning$value)))

        data_for_cleaning[, value := as.numeric(value)]
        n_clean_obs <- data_for_cleaning[!is.na(value), .N]
        n_clean_ppl <- data_for_cleaning[!is.na(value), .N, id][, .N]
        n_miss_num_conv_ppl <- n_total_ppl - n_clean_ppl

        # distribution
        pre_cutoff_sum <- data_for_cleaning[, .(
            pre_min_cutoff = quantile(value, 0.25, na.rm = TRUE) - 1.5 * IQR(value, na.rm = TRUE),
            pre_first_quartile = quantile(value, 0.25, na.rm = TRUE),
            pre_median = median(value, na.rm = TRUE),
            pre_mean = mean(value, na.rm = TRUE),
            pre_third_quartile = quantile(value, 0.75, na.rm = TRUE),
            pre_max_cutoff = quantile(value, 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE)
        )]

        # counts after cutoff
        n_miss_extreme_obs <- sum(is.na(exprs_clean(data_for_cleaning$value)))
        data_for_cleaning[, value := exprs_clean(value)]
        n_final_obs <- data_for_cleaning[!is.na(value), .N]
        n_final_ppl <- data_for_cleaning[!is.na(value), .N, id][, .N]
        n_miss_extreme_ppl <- n_clean_ppl - n_final_ppl
    })
    # summary stats
    final_sum <- data_for_cleaning[, .(
        final_min = min(value, na.rm = TRUE),
        final_first_quartile = quantile(value, 0.25, na.rm = TRUE),
        final_median = median(value, na.rm = TRUE),
        final_mean = mean(value, na.rm = TRUE),
        final_third_quartile = quantile(value, 0.75, na.rm = TRUE),
        final_max = max(value, na.rm = TRUE)
    )]

    # median and mean
    med_name <- paste0(name, "_med")
    mn_name <- paste0(name, "_mn")
    n_name <- paste0(name, "_n")
    data_for_cleaning_cleaned <- data_for_cleaning[!is.na(value), .(
        data_for_cleaning_med = median(value, na.rm = TRUE),
        data_for_cleaning_mn = mean(value, na.rm = TRUE),
        data_for_cleaning_obs = .N
    ), by = id] |>
        dplyr::rename(
            {{ med_name }} := data_for_cleaning_med,
            {{ mn_name }} := data_for_cleaning_mn,
            {{ n_name }} := data_for_cleaning_obs
        )

    cleaning_summary <- cbind(
        data.table(
            loinc_code = loinc_code,
            loinc_name = name,
            n_total_obs = n_total_obs,
            n_total_ppl = n_total_ppl,
            n_miss_num_conv_obs = n_miss_num_conv_obs,
            n_clean_obs = n_clean_obs,
            n_clean_ppl = n_clean_ppl,
            n_miss_num_conv_ppl = n_miss_num_conv_ppl
        ),
        pre_cutoff_sum,
        data.table(
            n_miss_extreme_obs = n_miss_extreme_obs,
            n_final_obs = n_final_obs,
            n_final_ppl = n_final_ppl,
            n_miss_extreme_ppl = n_miss_extreme_ppl,
            n_loss = n_total_ppl - n_final_ppl
        ),
        final_sum
    )

    return(
        list(
            data = data_for_cleaning_cleaned,
            cleaning_summary = cleaning_summary
        )
    )
}

