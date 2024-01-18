# Script for curating demographic, social hx, anthropometric, and lab data
# max salvatore
# 2024-01-15

# LIBRARIES --------------------------------------------------------------------
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

source("fn/curate_data_fn.R")

# DEMOGRAPHICS DATA ------------------------------------------------------------
cli_progress_step("Demographics")
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

keep_these_demo_vars <- c("id", "age", "age_enroll", "sex", "race_eth", "married", "female")

demo <- unique(demo[, ..keep_these_demo_vars])
qsave(
    x        = demo,
    file     = glue("data/private/mgi_demo_{version}.qs"),
    nthreads = 4
)

# SOCIAL HISTORY ---------------------------------------------------------------
cli_progress_step("Smoking and alcohol")
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

qsave(
    x        = smoking,
    file     = glue("data/private/mgi_smoking_{version}.qs"),
    nthreads = 4
)
qsave(
    x        = alcohol,
    file     = glue("data/private/mgi_alcohol_{version}.qs"),
    nthreads = 4
)
qsave(
    x        = social,
    file     = glue("data/private/mgi_social_{version}.qs"),
    nthreads = 4
)

rm(social)

# ANTHROPOMETRICS --------------------------------------------------------------
cli_progress_step("Anthropometrics")
anthro <- fread(paste0(path, "Anthropometrics_2023-03-22.txt")) |>
    clean_names()
setnames(
    anthro,
    old = c("de_id_patient_id"),
    new = c("id")
)
anthro <- anthro[, .(id, weight_kg, height_cm, bmi, dsb = days_since_birth)][, ysb := round(dsb / 365.25, 1)]

anthro[, `:=` (
    weight_kg_sd4 = sd4_clean(weight_kg),
    height_cm_sd4 = sd4_clean(height_cm),
    bmi_sd4 = sd4_clean(bmi),
    weight_kg_exprs = exprs_clean(weight_kg),
    height_cm_exprs = exprs_clean(height_cm),
    bmi_exprs = exprs_clean(bmi)
)]

ht_wt <- anthro[, .(
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

### thinness
ht_wt[, thinness_exprs_med := fifelse(
    bmi_exprs_med <= 19, 1, 0
)]

###

qsave(
    x        = anthro,
    file     = glue("data/private/mgi_anthro_{version}.qs"),
    nthreads = 4
)
qsave(
    x        = ht_wt,
    file     = glue("data/private/mgi_ht_wt_{version}.qs"),
    nthreads = 4
)

rm(anthro)

# VITALS -----------------------------------------------------------------------
cli_progress_step("Vitals")
vitals <- fread(paste0(path, "EncounterVitals_", as.Date(version, "%Y%m%d"), ".txt"),
    fill = TRUE, colClasses = "character"
) |>
    clean_names()
setnames(vitals, c("days_since_birth", "de_id_patient_id"), c("dsb", "id"))
vitals[, value := stringi::stri_trans_general(value, "latin-ascii")]

sbp <- vitals[term_subset == "BP Systolic", ][, value := as.numeric(value)][!is.na(value), ]
dbp <- vitals[term_subset == "BP Diastolic", ][, value := as.numeric(value)][!is.na(value), ]

sbp <- sbp[, value := exprs_clean(value)][!is.na(value), ]
dbp <- dbp[, value := exprs_clean(value)][!is.na(value), ]

sbp <- sbp[!is.na(value), .(
    sbp_med = median(value, na.rm = TRUE),
    sbp_mn = mean(value, na.rm = TRUE),
    sbp_obs = .N
), by = id]

dbp <- dbp[!is.na(value), .(
    dbp_med = median(value, na.rm = TRUE),
    dbp_mn = mean(value, na.rm = TRUE),
    dbp_obs = .N
), by = id]

bp <- merge.data.table(sbp, dbp, by = "id", all = TRUE)[, in_vitals := 1]

qsave(
    x        = bp,
    file     = glue("data/private/mgi_vitals_{version}.qs"),
    nthreads = 4
)

## combine
#
ht_wt_keep <- c("id", grep("_exprs_med", names(ht_wt), value = TRUE))
ht_wt <- ht_wt[, ..ht_wt_keep]
names(ht_wt) <- gsub("_exprs", "", names(ht_wt))

bp_keep <- c("id", grep("_med", names(bp), value = TRUE))
bp <- bp[, ..bp_keep]

demo_combined <- reduce(
    list(demo, smoking, alcohol, ht_wt, bp),
    full_join,
    by = "id"
)
keep_demo_combined <- names(demo_combined)[!grepl("_obs", names(demo_combined))]
demo_combined <- demo_combined[, ..keep_demo_combined]

qsave(
    x        = demo_combined,
    file     = glue("data/private/mgi_demo_comb_{version}.qs"),
    nthreads = 4
)

# LABS -------------------------------------------------------------------------
cli_progress_step("Labs")
# lab results file path
labs_results_file <- paste0(path, "LabResults_", as.Date(version, "%Y%m%d"), ".txt")

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

### total cholesterol
hdl <- labs[loinc == "2085-9", ][, value := as.numeric(value)]
ldl <- labs[loinc == "13457-7", ][, value := as.numeric(value)]
tri <- labs[loinc == "2571-8", ][, value := as.numeric(value)]

hdl[, .N, by = c("id", "collection_dsb")][N > 1, ]
ldl[, .N, by = c("id", "collection_dsb")][N > 1, ]
tri[, .N, by = c("id", "collection_dsb")][N > 1, ]

remove_any_duplicate <- function(x, cols) {
    x[!(duplicated(x[, ..cols]) | duplicated(x[, ..cols], fromLast = TRUE)), ][]
}

hdl <- remove_any_duplicate(hdl, c("id", "collection_dsb"))
ldl <- remove_any_duplicate(ldl, c("id", "collection_dsb"))
tri <- remove_any_duplicate(tri, c("id", "collection_dsb"))

chol <- reduce(
    list(
        hdl[, .(id, dsb = collection_dsb, hdl = value)],
        ldl[, .(id, dsb = collection_dsb, ldl = value)],
        tri[, .(id, dsb = collection_dsb, tri = value)]
    ),
    merge.data.table,
    by = c("id", "dsb"),
    all = TRUE
) |>
    na.omit()

chol[, `:=`(
    total_chol = hdl + ldl + (0.2 * tri)
)]

chol[, `:=`(
    total_chol = exprs_clean(total_chol)
)]

total_chol <- chol[!is.na(total_chol), .(
    total_chol_med = median(total_chol, na.rm = TRUE),
    total_chol_mn = mean(total_chol, na.rm = TRUE),
    total_chol_obs = .N
), by = id]

cleaned_labs <- merge.data.table(
    cleaned_labs,
    total_chol,
    by = "id",
    all = TRUE
)
###

qsave(
    x        = cleaned_labs,
    file     = glue("data/private/mgi_labs_cleaned_{version}.qs"),
    nthreads = 4
)

medians_only <- c("id", grep("_med", names(cleaned_labs), value = TRUE))
cleaned_labs_med <- cleaned_labs[, ..medians_only]
qsave(
    x        = cleaned_labs_med,
    file     = glue("data/private/mgi_labs_cleaned_medians_{version}.qs"),
    nthreads = 4
)

# summaries
fwrite(
    x = rbindlist(map(labs_cleaned, 2)),
    file = glue("data/private/mgi_labs_cleaning_summary_{version}.csv")
)

##
demo_labs <- full_join(
    demo_combined,
    cleaned_labs_med,
    by = "id"
)

qsave(
    x        = demo_labs,
    file     = glue("data/private/mgi_demo_comb_labs_{version}.qs"),
    nthreads = 4
)

cli::cli_alert_success("Done! 🎉")
