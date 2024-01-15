ms::libri(ms, data.table, tidyverse, mice, cli, glue, janitor, qs)

version <- "20230322"
path <- paste0("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/", version, "/")
# Freeze 3

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

# # read in LOINC codes to filter
# loinc_code_file <- read_delim("data/public/loinc_codes.csv", show_col_types = FALSE) |>
#     clean_names()
# loinc_codes <- loinc_code_file |>
#     dplyr::filter(ex_prs == 1) |>
#     dplyr::pull(loinc_code) |>
#     unique()

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
#             dplyr::filter(loinc %in% loinc_codes)
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

labs <- left_join(
    labs,
    loinc_code_file |>
        dplyr::filter(ex_prs == 1) |>
        dplyr::select(loinc = loinc_code, test_order_name, exposure),
    by = "loinc"
)


hdl <- labs[loinc == "2085-9", ]
n_total_obs <- hdl[, .N]
n_total_ppl <- hdl[, .N, id][, .N]
n_miss_num_conv_obs <- sum(is.na(as.numeric(hdl$value)))

hdl[, value := as.numeric(value)]
n_clean_obs <- hdl[!is.na(value), .N]
n_clean_ppl <- hdl[!is.na(value), .N, id][, .N]
n_miss_num_conv_ppl <- n_total_ppl - n_clean_ppl

pre_cutoff_sum <- hdl[, .(
    min_cutoff = quantile(value, 0.25, na.rm = TRUE) - 1.5 * IQR(value, na.rm = TRUE),
    first_quartile = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    third_quartile = quantile(value, 0.75, na.rm = TRUE),
    max_cutoff = quantile(value, 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE)
)]


n_miss_extreme_obs <- sum(is.na(exprs_clean(hdl$value)))
hdl[, value := exprs_clean(value)]
n_final_obs <- hdl[!is.na(value), .N]
n_final_ppl <- hdl[!is.na(value), .N, id][, .N]
n_miss_extreme_ppl <- n_clean_ppl - n_final_ppl

final_sum <- hdl[, .(
    min = min(value, na.rm = TRUE),
    first_quartile = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    third_quartile = quantile(value, 0.75, na.rm = TRUE),
    max = max(value, na.rm = TRUE)
)]

hdl_cleaned <- hdl[!is.na(value), .(
    hdl_med = median(value, na.rm = TRUE),
    hdl_mn = mean(value, na.rm = TRUE),
    hdl_obs = .N
), by = id]

hdl 

length(unique(hdl$id))
summary(hdl[, .N, id][, N])
summary(hdl[, value])



labs |>
    group_by(loinc) |>
    (\(x) x[is.na(as.numeric(value)), ])() |>
    (\(x) x[, unique(value)])()


sapply(labs, \(x) length(unique(x)))