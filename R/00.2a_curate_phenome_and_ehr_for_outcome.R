# libraries, paths, and such ---------------------------------------------------
ms::libri(
  data.table, glue, qs, parallel, optparse, PheWAS/PheWAS, cli, ms, comorbidity,
  labelled
)

# optparse list ----
option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

cli_alert(glue("using mgi cohort version {opt$mgi_version}"))

out_path <- glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_three/data/private/mgi/{opt$mgi_version}/")
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

# source("fn/files-utils.R")
# file_paths <- get_files(mgi_version = opt$mgi_version)

# load data --------------------------------------------------------------------
cli_alert("loading data...")
# load demo
# demo <- qread(glue("data/private/mgi_demo_comb_labs_20230322.qs"))

# ### MAP ICD9/ICD10 codes
# mgi_icd9  <- get(load("/net/junglebook/magic_data/EHRdata/20220822/phenomes/UNFILTERED_20220822/UNFILTERED_20220822_ICD9_Phecodes_Birthyears.Rsav"))
# mgi_icd10 <- get(load("/net/junglebook/magic_data/EHRdata/20220822/phenomes/UNFILTERED_20220822/UNFILTERED_20220822_ICD10_Phecodes_Birthyears.Rsav"))

# mgi_icd <- dplyr::bind_rows(
#   mgi_icd9 |>
#     dplyr::mutate(vocabulary_id = "ICD9CM") |>
#     dplyr::select(IID, vocabulary_id, DiagnosisCode, DaysSinceBirth),
#   mgi_icd10 |>
#     dplyr::mutate(vocabulary_id = "ICD10CM") |>
#     dplyr::select(IID, vocabulary_id, DiagnosisCode, DaysSinceBirth)
# )

# in_demo_and_phenome <- intersect(unique(demo[, id]), unique(mgi_icd[, IID]))

# rm(mgi_icd9, mgi_icd10)

# qsave(
#   mgi_icd |> 
#     dplyr::select(IID, lexicon = vocabulary_id, DiagnosisCode, DaysSinceBirth),
#   file = glue("{out_path}MGI_ICD_{opt$mgi_version}.qs")
# )

# mgi_icd <- mgi_icd |>
#   dplyr::rename(id = IID, code = DiagnosisCode, dsb = DaysSinceBirth)

# ## external
# ### mapping table
# phecodex_map <- ms::phecodex_icdcm_map
# ### rollup map
# phecodex_rollup <- ms::phecodex_rollup
# ### sex restriction map
# phecodex_sex <- ms::phecodex_sex
# ### phecodex info
# phecodex_info <- ms::phecodex_labels

# # map phenotypes ---------------------------------------------------------------
# mapped <- mapCodesToPhecodes(
#   input = mgi_icd, vocabulary.map = phecodex_map,
#   rollup.map = phecodex_rollup, make.distinct = TRUE
# )

# id_sex <- demo |> 
#   dplyr::select(id, sex)

# # remove sex-specific phecodes that are discordant with individual's reported sex at birth
# restrictPhecodesBySex_mod <- function(phenotypes, id.sex, by_var = "person_id", sex_var = "sex") {
#   data <- merge.data.table(
#     as.data.table(phenotypes),
#     as.data.table(id.sex),
#     by = by_var, all.x = TRUE
#   )
#   # Get the restrictions found in the phenotypes data frame
#   current_sex_restriction <- phecodex_sex[phecodex_sex$phecode %in% unique(data[, phecode]), ] |>
#     as.data.table()
#   # Get male and female-only phenotypes
#   male_only <- current_sex_restriction[current_sex_restriction$male_only, phecode]
#   female_only <- current_sex_restriction[current_sex_restriction$female_only, phecode]
#   # Set row column matches to NA where inds of a sex meet restricted phenotypes
#   data[phecode %in% male_only & sex == "Female", phecode := NA]
#   data[phecode %in% female_only & sex == "Male", phecode := NA]

#   na.omit(data)[, (sex_var) := NULL][]
# }

# ## remove sex person-phenotype discordant pairs
# restricted <- restrictPhecodesBySex_mod(mapped, id.sex = id_sex, by_var = "id")

# in_demo_and_phenome <- intersect(unique(demo[, id]), unique(restricted[, id]))

# restricted <- restricted[id %in% in_demo_and_phenome, ]

## save
restricted <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs"))

# get first phecode per person -------------------------------------------------
first_restricted <- restricted[restricted[, .I[which.min(dsb)], by = c("id", "phecode")][["V1"]], ]

mgi_icd <- qread(glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_three/data/private/mgi/{opt$mgi_version}/MGI_ICD_{opt$mgi_version}.qs"))

# outcome ----------------------------------------------------------------------
outcome <- "EM_202.2"
cases <- first_restricted[phecode == outcome, unique(id)]
case_restricted    <- first_restricted[id %in% cases, ]
control_restricted <- first_restricted[!(id %in% cases), ]

case_index      <- case_restricted[phecode == outcome, .(id, index_dsb = dsb)]

qsave(
  case_index,
  file = glue("data/private/mgi/{opt$mgi_version}/MGI_{outcome}_INDEX_{opt$mgi_version}.qs")
)

case_restricted <- merge.data.table(case_restricted, case_index, by = "id", all.x = TRUE)
case_sub <- case_restricted[dsb < index_dsb, ]

stack_phe <- bind_rows(case_sub, control_restricted) |> as.data.table()
stack_ids <- unique(stack_phe[, id])

# subset ICD code data for CCI calculation later
icd_sub <- mgi_icd[, .(id = IID, code = DiagnosisCode, dsb = DaysSinceBirth, vocabulary_id = lexicon)]
icd_sub <- icd_sub[id %in% stack_phe[, unique(id)], ]

# pim --------------------------------------------------------------------------
pim <- dcast(
  stack_phe[, .(id, phecode = phecode)],
  id ~ phecode,
  value.var = "phecode",
  fun.aggregate = length,
  fill = 0
)
pim[, case := fifelse(id %in% cases, 1, 0)]
# collect integer variable names from data.table
int_vars <- names(pim)[which(sapply(pim, is.integer))]
pim[, (int_vars) := lapply(.SD, \(x) as.numeric(x > 0)), .SDcols = int_vars]
int_vars_sub <- int_vars[int_vars != "case"]
for (i in seq_along(int_vars_sub)) {
  var_label(pim[[int_vars_sub[i]]]) <- ms::pheinfox[phecode == int_vars_sub[i], description]
}

qsave(
  pim,
  file = glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}_{outcome}.qs")
)

### first and last dsb
restricted_sub <- restricted[id %in% stack_ids, ]
first_dsb <- unique(restricted_sub[restricted_sub[
  ,
  .I[which.min(dsb)],
  "id"
][["V1"]]][
  ,
  .(
    id,
    first_dsbx = dsb,
    age_at_first_diagnosisx = round(dsb / 365.25, 1)
  )
])
last_dsb <- unique(restricted_sub[restricted_sub[, .I[which.max(dsb)], "id"][["V1"]]][, .(
  id,
  last_dsbx = dsb,
  age_at_last_diagnosisx = round(dsb / 365.25, 1)
)])

demo <- merge.data.table(
  first_dsb,
  last_dsb, by = "id", all = TRUE
)[, length_followupx := round((last_dsbx - first_dsbx) / 365.25)]

###
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")[
  keep == 1, phecode
]

comorbid <- list(
  "cadx"                = list("phecodes" = c("CV_404.2")),
  "diabetesx"           = list("phecodes" = c("EM_202")),
  "hypertensionx"       = list("phecodes" = c("CV_401")),
  "mixed_hypertensionx" = list("phecodes" = c("EM_239.3")),
  "vitamin_dx"          = list("phecodes" = c("EM_232.4")),
  "depressionx"         = list("phecodes" = c("MB_286.2")),
  "anxietyx"            = list("phecodes" = c("MB_288")),
  "bipolarx"            = list("phecodes" = c("MB_286.1")),
  "cancerx"             = list("phecodes" = cancer_phecodesx)
)

# identify cases and create indicator variables --------------------------------
cli_alert("identifying cases and creating indicator variables...")
cli_progress_bar("deriving comorbidities", total = length(names(comorbid)))
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- restricted_sub[phecode %in% comorbid[[i]][["phecodes"]], id] |>
    unique()
  cli_progress_update()
}
cli_progress_done()

cli_progress_bar("deriving comorbidity indicator variables", total = length(names(comorbid)))
for (i in names(comorbid)) {
  set(demo, j = i, value = fifelse(demo[["id"]] %in% comorbid[[i]][["ids"]], 1, 0))
  cli_progress_update()
}
cli_progress_done()

demo[, triglyceridesx := fifelse(hypertensionx == 0 & mixed_hypertensionx == 0, 0, 1)]

demo[, age_catx := between(age_at_last_diagnosisx, 0, 5.99) +
  2 * between(age_at_last_diagnosisx, 6, 11.99) +
  3 * between(age_at_last_diagnosisx, 12, 19.99) +
  4 * between(age_at_last_diagnosisx, 20, 39.99) +
  5 * between(age_at_last_diagnosisx, 40, 59.99) +
  6 * between(age_at_last_diagnosisx, 60, 150.99)]

###
# ADD INDICATOR FOR FIRST CANCER DIAGNOSIS
cancer_phecodesx_spec <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")[
  keep == 1 & specific == 1, phecode
]
first_cancer <- stack_phe[phecode %in% cancer_phecodesx_spec, ]
first_cancer <- first_cancer[
  first_cancer[, .I[dsb == min(dsb, na.rm = TRUE)], id][["V1"]]
]

mgi_first_cancer <- data.table(
  id = first_cancer[, unique(id)]
)

first_cancer_ids <- list()
for (i in seq_along(cancer_phecodesx_spec)) {
  first_cancer_ids[[i]] <- first_cancer[phecode == cancer_phecodesx_spec[i], unique(id)]
  names(first_cancer_ids)[i] <- cancer_phecodesx_spec[i]
}

for (i in seq_along(cancer_phecodesx_spec)) {
  set(
    mgi_first_cancer,
    j     = paste0(cancer_phecodesx_spec[i], "_first"),
    value = fifelse(mgi_first_cancer[["id"]] %in% first_cancer_ids[[i]], 1, 0)
  )
}

### comorbidity scores
# charlson comorbidity index
charlson10 <- comorbidity(icd_sub[vocabulary_id == "ICD10CM", ], id = "id", code = "code", map = "charlson_icd10_quan", assign0 = FALSE)
charlson9  <- comorbidity(icd_sub[vocabulary_id == "ICD9CM", ], id = "id", code = "code", map = "charlson_icd9_quan", assign0 = FALSE)

charlson <- merge.data.table(
  charlson10,
  charlson9,
  by = "id",
  all = TRUE
) |> as.data.table()

charl_comorbid_names <- gsub(".x", "", grep(".x", names(charlson), value = TRUE))
for (i in seq_along(charl_comorbid_names)) {
  charlson[, (charl_comorbid_names[i]) := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c(paste0(charl_comorbid_names[i], ".x"), paste0(charl_comorbid_names[i], ".y"))]
}
keep_charl                        <- c("id", charl_comorbid_names)
charlson                          <- charlson[, ..keep_charl]
attr(charlson, "class")           <- attributes(charlson10)$class
attr(charlson, "variable.labels") <- attributes(charlson10)$variable.labels
attr(charlson, "map")             <- attributes(charlson10)$map
un_charl_scores                   <- score(charlson, weights = NULL, assign0 = FALSE)
we_charl_scores                   <- score(charlson, weights = "quan", assign0 = FALSE)

charl_scores <- data.table(
  id = charlson$id,
  charlson_unweighted = un_charl_scores,
  charlson_quan_weighted = we_charl_scores
)

# elixhauser comorbidity index
elixhauser10 <- comorbidity(icd_sub[vocabulary_id == "ICD10CM", ], id = "id", code = "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
elixhauser9  <- comorbidity(icd_sub[vocabulary_id == "ICD9CM", ], id = "id", code = "code", map = "elixhauser_icd9_quan", assign0 = FALSE)

elixhauser <- merge.data.table(
  elixhauser10,
  elixhauser9,
  by = "id",
  all = TRUE
) |> as.data.table()

elix_comorbid_names <- gsub(".x", "", grep(".x", names(elixhauser), value = TRUE))
for (i in seq_along(elix_comorbid_names)) {
  elixhauser[, (elix_comorbid_names[i]) := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c(paste0(elix_comorbid_names[i], ".x"), paste0(elix_comorbid_names[i], ".y"))]
}
keep_elix                           <- c("id", elix_comorbid_names)
elixhauser                          <- elixhauser[, ..keep_elix]
attr(elixhauser, "class")           <- attributes(elixhauser10)$class
attr(elixhauser, "variable.labels") <- attributes(elixhauser10)$variable.labels
attr(elixhauser, "map")             <- attributes(elixhauser10)$map
un_elix_scores                      <- score(elixhauser, weights = NULL, assign0 = FALSE)
we_elix_scores                      <- score(elixhauser, weights = "vw", assign0 = FALSE)

elix_scores <- data.table(
  id = elixhauser$id,
  elixhauser_unweighted = un_elix_scores,
  elixhauser_vw_weighted = we_elix_scores
)

comorbid_scores <- merge.data.table(
  charl_scores,
  elix_scores,
  by = "id",
  all = TRUE
) |> as.data.table()

demo <- merge.data.table(
  demo,
  comorbid_scores,
  by = "id",
  all = TRUE
)
###


###
qsave(
  as_tibble(charlson),
  file = glue("data/private/mgi/{opt$mgi_version}/charlsonx_{opt$mgi_version}_{outcome}.qs")
)
qsave(
  as_tibble(elixhauser),
  file = glue("data/private/mgi/{opt$mgi_version}/elixhauserx_{opt$mgi_version}_{outcome}.qs")
)

qsave(
  mgi_first_cancer,
  file = glue("data/private/mgi/{opt$mgi_version}/first_cancerx_{opt$mgi_version}_{outcome}.qs")
)
###

qsave(
  demo,
  file = glue("data/private/mgi/{opt$mgi_version}/ehr_{opt$mgi_version}_{outcome}.qs")
)

cli_alert_success("Finished! 🎉")
