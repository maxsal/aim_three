# Script for identifying individuals in all datasets
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

# demo
demo <- qread(glue("data/private/mgi_demo_comb_labs_{version}.qs"))[, in_demo := 1]
demo_vars <- c("age", "age_enroll", "sex", "race_eth", "married", "female")
smoke_vars <- c("smoke_ever")
alcohol_vars <- c("alcohol_ever")
ht_wt_vars <- c("weight_kg_med", "heigh_cm_med", "bmi_med")

loinc_code_file <- fread("data/public/loinc_codes.csv") |>
    clean_names() |>
    dplyr::filter(ex_prs == 1) |>
    dplyr::mutate(
        exposure = make_clean_names(exposure)
    )

labs_vars <- names(demo)[names(demo) %in% paste0(loinc_code_file[, exposure], "_med")]

demo[apply(demo[, ..labs_vars], MARGIN = 1, \(x) any(!is.na(x))), in_labs := 1]
demo <- demo[!is.na(sex), ]

# diag
diag <- qread(glue("data/private/mgi_pim0_{version}.qs"))[, in_diag := 1]

# geno
geno <- fread(glue("data/private/mgi_geno_ids_{version}.txt"))[, in_geno := 1]

keep_ids <- reduce(
    list(demo[complete.cases(demo[, .(age, sex)]), id], diag[, id], geno[, id]),
    intersect
)

comb <- reduce(
    list(demo, diag, geno),
    merge.data.table, by = "id", all = TRUE
)

comb[, .N, in_demo]
comb[ is.na(in_demo), ]

# restrict to individuals non-missing age and sex and are present in phenome and genome
demog <- demo[id %in% keep_ids, ]
diagg <- diag[id %in% keep_ids, ]
int_vars <- names(diagg)[sapply(diagg, is.numeric)]
diagg[, (int_vars) := lapply(.SD, \(x) {
    ifelse(
        is.na(x),
        0,
        as.numeric(x > 0)
    )
}), .SDcols = int_vars]
genog <- geno[id %in% keep_ids, ]

combined <- reduce(
    list(demog, diagg, genog),
    merge.data.table, by = "id", all = TRUE
)

n <- combined[, length(unique(id))]
not_phenome_vars <- names(combined)[!(names(combined) %in% names(diag))]
miss_count <- sapply(combined[, ..not_phenome_vars], function(x) sum(is.na(x)))
missing_summary <- data.table(
    var = names(miss_count),
    n_miss = miss_count,
    pct_miss = round(miss_count / n * 100, 1)
)[order(-pct_miss), ]

fwrite(
    missing_summary,
    glue("data/private/missing_summary_exprs_traits_{version}.csv")
)

demog[complete.cases(demog), ]