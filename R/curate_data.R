# Script for curating demographic, social hx, anthropometric, and lab data
# max salvatore
# 2024-01-15
ms::libri(ms, data.table, tidyverse, mice, cli, glue, janitor, qs, optparse)

option_list <- list(
    make_option(c("-v", "--version"),
        type = "character", default = "20230322",
        help = "MGI MAGIC data pull version [default = %default]"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

version <- opt$version
path <- paste0("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/", version, "/")

# data
## demographics
demo <- fread(paste0(path, "Demographics_", as.Date(version, "%Y%m%d"), ".txt")) |>
    clean_names()
setnames(
    demo,
    old = c("de_id_patient_id"),
    new = c("id")
)
demo[, `:=` (
    sex = fcase(
        gender_code == "M", "Male",
        gender_code == "F", "Female"
    ),
    married = fcase(
        marital_status_code == "M", "Married",
        marital_status_code == "U", "Unmarried"
    ),
    race_eth = fcase(
        ethnicity_name == "Non-Hispanic or Latino" & race_name == "Caucasian", "NHW",
        ethnicity_name == "Non-Hispanic or Latino" & race_name == "African American", "NHB",
        ethnicity_name == "Non-Hispanic or Latino" & race_name == "Asian", "NHA",
        ethnicity_name == "Hispanic or Latino", "Hispanic",
        race_name == "Other", "Other"
    ),
    age_enroll = round(enrollment_days_since_birth / 365.25, 1)
)]
demo[, female := fcase(
    sex == "Female", 1,
    sex == "Male", 0
)]

keep_these_demo_vars <- c("id", "age", "age_enroll", "sex", "race_eth", "married")

## SocialHx
social <- fread(paste0(path, "SocialHx_", as.Date(version, "%Y%m%d"), ".txt")) |>
    clean_names()
setnames(
    social,
    old = c("de_id_patient_id"),
    new = c("id")
)

smoking <- social[, .(id, smoking_status)] |>
    dplyr::summarize(
        smoke_ever = fcase(
            any(smoking_status %in% c("Current", "Former")), 1,
            any(smoking_status %in% c("Never")), 0
        ),
        smoke_total_obs = n(),
        smoke_ans_obs   = sum(!is.na(smoking_status)),
        .by = id
    ) |> as.data.table()

alcohol <- social[, .(id, alcohol_use_status)] |>
    dplyr::summarize(
        alcohol_ever = fcase(
            any(alcohol_use_status %in% c("Yes")), 1,
            any(alcohol_use_status %in% c("No")), 0
        ),
        alcohol_total_obs = n(),
        alcohol_ans_obs   = sum(!is.na(alcohol_use_status)),
        .by = id
    ) |> as.data.table()

## Anthropometrics
anthro <- fread(paste0(path, "Anthropometrics_2023-03-22.txt")) |>
    clean_names()
setnames(
    anthro,
    old = c("de_id_patient_id"),
    new = c("id")
)
anthro <- anthro[, .(id, weight_kg, height_cm, bmi, dsb = days_since_birth)][, ysb := round(dsb / 365.25, 1)]

# approach from Li, Chen & Moore 2019 (PMID: 31329892)
sd4_clean <- function(x) {
    mn <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    x[x > mn + 4 * sd] <- NA
    x[x < mn - 4 * sd] <- NA
    x
}
# approach from Ma,...,Fritsche 2022 
exprs_clean <- function(x) {
    x_qt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
    x_iqr <- IQR(x, na.rm = TRUE)
    x[(x > (x_qt[2] + 1.5*x_iqr))] <- NA
    x[(x < (x_qt[1] - 1.5*x_iqr))] <- NA
    x
}

anthro[, `:=` (
    weight_kg_sd4 = sd4_clean(weight_kg),
    height_cm_sd4 = sd4_clean(height_cm),
    bmi_sd4 = sd4_clean(bmi),
    weight_kg_exprs = exprs_clean(weight_kg),
    height_cm_exprs = exprs_clean(height_cm),
    bmi_exprs = exprs_clean(bmi)
)]

ht_wt <- tmp[, .(
    weight_kg_sd4_med = median(weight_kg_sd4, na.rm = TRUE),
    weight_kg_sd4_mn  = mean(weight_kg_sd4, na.rm = TRUE),
    height_cm_sd4_med = median(height_cm_sd4, na.rm = TRUE),
    height_cm_sd4_mn  = mean(height_cm_sd4, na.rm = TRUE),
    bmi_sd4_med = median(bmi_sd4, na.rm = TRUE),
    bmi_sd4_mn  = mean(bmi_sd4, na.rm = TRUE),

    weight_kg_exprs_med = median(weight_kg_exprs, na.rm = TRUE),
    weight_kg_exprs_mn = mean(weight_kg_exprs, na.rm = TRUE),
    height_cm_exprs_med = median(height_cm_exprs, na.rm = TRUE),
    height_cm_exprs_mn = mean(height_cm_exprs, na.rm = TRUE),
    bmi_exprs_med = median(bmi_exprs, na.rm = TRUE),
    bmi_exprs_mn = mean(bmi_exprs, na.rm = TRUE),
    bmi_n = sum(!is.na(bmi))
), id]

## Labs
# # lab results file path
# labs_results_file <- paste0(path, "LabResults_", as.Date(version, "%Y%m%d"), ".txt")

# # helper function to read number of rows in a file
# get_line_count <- function(file) {
#     f <- file(file, open = "rb")
#     nlines <- 0L
#     while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
#         nlines <- nlines + sum(chunk == as.raw(10L))
#     }
#     nlines
# }

# read in LOINC codes to filter
loinc_code_file <- read_delim("data/public/loinc_codes.csv", show_col_types = FALSE) |>
    clean_names() |> as.data.table()
loinc_codes <- loinc_code_file |>
    dplyr::filter(ex_prs == 1)
loinc_codes_only <- loinc_codes[, unique(loinc_code)]

# # prepare for reading in chunks
# lines_per_chunk <- 1e6 # 1 million lines per chunk
# total_lines     <- get_line_count(labs_results_file) # ~72 million lines
# chunks          <- ceiling(total_lines / lines_per_chunk) # ~72 chunks

# var_names <- read_delim(labs_results_file, n_max = 1, show_col_types = FALSE) |>
#     clean_names() |>
#     names()

# # read in lab results file in chunks of 1 million lines
# lab_res <- map(
#     1:chunks,
#     \(chunk) {
#             read_delim(labs_results_file,
#                 skip = 1 + (chunk - 1) * lines_per_chunk, n_max = lines_per_chunk,
#                 col_names = var_names, col_types = cols(.default = "c"),
#                 progress = FALSE
#             ) |>
#             dplyr::filter(loinc %in% loinc_codes_only)
#     },
#     .progress = TRUE
# )

# # combine chunks into one data.table
# labs <- bind_rows(lab_res) |>
#     select(
#         id = de_id_patient_id,
#         collection_dsb = collection_date_days_since_birth,
#         order_code, order_name,
#         result_code, result_name,
#         loinc, value, unit,
#         range, specimen_source, hilonormal_flag
#     ) |>
#     as.data.table()


# labs[, length(unique(id))]

# qsave(
#     x         = labs,
#     file      = glue("data/private/mgi_labs_raw_{version}.qs"),
#     nthreads  = 4
# )

labs <- qread(glue("data/private/mgi_labs_raw_{version}.qs"))

labs <- left_join(
    labs,
    loinc_code_file |>
        dplyr::filter(ex_prs == 1) |>
        dplyr::select(loinc = loinc_code, test_order_name, exposure),
    by = "loinc"
)

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
            {{med_name}} := data_for_cleaning_med,
            {{mn_name}} := data_for_cleaning_mn,
            {{n_name}} := data_for_cleaning_obs
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

clean_these_labs        <- loinc_codes[, loinc_code]
names(clean_these_labs) <- janitor::make_clean_names(loinc_codes[, exposure])

labs_cleaned <- map(
    seq_along(clean_these_labs),
    \(i) {
        clean_labs(
            lab_data   = labs,
            loinc_code = clean_these_labs[i],
            name       = names(clean_these_labs)[i]
        )
    },
    .progress = TRUE
) |>
    set_names(names(clean_these_labs))

cleaned_labs <- reduce(map(labs_cleaned, 1), full_join, by = "id")
qsave(
    x        = cleaned_labs,
    file     = glue("data/private/mgi_labs_cleaned_{version}.qs"),
    nthreads = 4
)

medians_only <- c("id", grep("_med", names(cleaned_labs), value = TRUE))
cleaned_labs_med <- cleaned_labs[, ..medians_only]
qsave(
    x        = cleaned_labs,
    file     = glue("data/private/mgi_labs_cleaned_medians_{version}.qs"),
    nthreads = 4
)

# summaries
fwrite(
    x = rbindlist(map(labs_cleaned, 2)),
    file = glue("data/private/mgi_labs_cleaning_summary_{version}.csv")
)

