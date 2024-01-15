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

# DIAGNOSTIC DATA --------------------------------------------------------------
diag_data_file <- paste0(path, "Diagnosis_", as.Date(version, "%Y%m%d"), ".txt")
diag_columns <- c("DeId_PatientID", "DaysSinceBirth", "DiagnosisCode", "Lexicon")
diag <- fread(diag_data_file, select = diag_columns, colClasses = "character")
setnames(diag, old = c("DeId_PatientID", "DaysSinceBirth", "DiagnosisCode", "Lexicon"),
    new = c("id", "dsb", "code", "vocabulary_id"))
diag[, dsb := as.numeric(dsb)]
diag[vocabulary_id == "ICD9", vocabulary_id := "ICD9CM"]
diag[vocabulary_id == "ICD10", vocabulary_id := "ICD10CM"]

## PHECODE 1.2 -----------------------------------------------------------------
phecodes <- PheWAS::mapCodesToPhecodes(
    input = diag
)

# read in demo 
demo <- qread(glue("data/private/mgi_demo_labs_{version}.qs"))[, .(
    id,
    sex = fcase(sex == "Male", "M", sex == "Female", "F")
)]

restrictPhecodesBySex_mod <- function(phenotypes, id.sex, id.var = "id", sex.var = "sex", sex.restriction = PheWAS::sex_restriction) {
    data <- merge.data.table(
        as.data.table(phenotypes),
        as.data.table(id.sex),
        by = id.var, all.x = TRUE)
    # Get the restrictions found in the phenotypes data frame
    data <- merge.data.table(
        data,
        sex.restriction,
        by = "phecode",
        all.x = TRUE
    )
    data <- data[!(sex == "F" & male_only == TRUE) & !(sex == "M" & female_only == TRUE), ]
    # Return everything, sans sex
    rm_cols <- c(sex.var, "male_only", "female_only")
    data[, (rm_cols) := NULL][]
}


phecodes_restricted <- restrictPhecodesBySex_mod(
    phenotypes = phecodes,
    id.sex  = demo
)

first_restricted <- phecodes_restricted[phecodes_restricted[, .I[which.min(dsb)], by = c("id", "phecode")][["V1"]]]

qsave(
    phecodes_restricted,
    glue("data/private/mgi_phe_dsb_{version}.qs"),
    nthreads = 4
)

qsave(
    first_restricted,
    glue("data/private/mgi_first_phe_dsb_{version}.qs"),
    nthreads = 4
)

# identify exclusions
sample_ids <- intersect(pim[, unique(id)], demo[, unique(id)])
females <- demo[sex == "F", unique(id)]
males   <- demo[sex == "M", unique(id)]

femaleTF <- sample_ids %in% females
maleTF <- sample_ids %in% males

inclusions <- first_restricted[, c("id", "phecode")]
inclusions <- split(inclusions$id, inclusions$phecode)

exclusions <- first_restricted[, c("id", "phecode")]
exclusions <- split(exclusions$id, exclusions$phecode)

pheinfo <- fread("https://raw.githubusercontent.com/umich-cphds/createUKBphenome/master/data/Phecode_Definitions_FullTable_Modified.txt",
    colClasses = "character")
phecode_exclude <- as.data.table(PheWAS::phecode_exclude)

phecodes <- pheinfo[, unique(phecode)]

# create two empty data.frames (one with exclusion criteria applied, one without apply exclusion criteria); add X to phecodes
phenoOut <- phenoOut0 <- data.table("id" = sample_ids, matrix(NA, ncol = length(phecodes), nrow = length(sample_ids)))
colnames(phenoOut) <- c("id", paste0("X", phecodes))
colnames(phenoOut0) <- c("id", paste0("X", phecodes))

pb <- txtProgressBar(min = 0, max = length(phecodes), initial = 0, style = 3)
for (i in seq_along(phecodes)) {
    phecode_remove <- ""
    .phecode <- phecodes[i]

    # collect phecodes to include from controls
    exclude_phecodes <- .phecode
    if (nrow(phecode_exclude[code == .phecode]) > 0) {
        exclude_phecodes <- c(exclude_phecodes, phecode_exclude[code == .phecode, exclusion_criteria])           
    }
    exclude_phecodes <- unique(exclude_phecodes[which(exclude_phecodes %in% pheinfo[, phecode])])

    # collected sample sets in inclusion and exclusion criteria
    pinclusions <- unlist(inclusions[.phecode])
    pexclusions <- unique(unlist(exclusions[exclude_phecodes]))

    # phecode based case control defitions (without sex filter)
    cases <- sample_ids %in% pinclusions
    controls <- !sample_ids %in% pexclusions

    # case control status
    ccstatus <- ccstatus0 <- rep(NA, nrow(phenoOut))
    if (pheinfo[phecode == .phecode, sex] == "Both") {
        ccstatus[controls] <- 0
        ccstatus[cases] <- ccstatus0[cases] <- 1
        ccstatus0[!cases] <- 0
    } else if (pheinfo[phecode == .phecode, sex == "Female"]) {
        ccstatus[controls & femaleTF] <- 0
        ccstatus[cases & femaleTF] <- ccstats0[cases & femaleTF] <- 1
        ccstatus0[!cases & femaleTF] <- 0
        checkFemales <- sample_ids[cases & !femaleTF]
        if(length(checkFemales) > 0) notFemale[[phecode]] <- data.table(id = checkFemales, phecode)
    } else if (pheinfo[phecode == .phecode, sex == "Male"]) {
        ccstatus[controls & maleTF] <- 0
        ccstatus[cases & maleTF] <- ccstats0[cases & maleTF] <- 1
        ccstatus0[!cases & maleTF] <- 0
        checkMales <- sample_ids[cases & !maleTF]
        if(length(checkMales) > 0) notMale[[phecode]] <- data.table(id = checkMales, phecode)
    }
    phenoOut[[paste0("X", .phecode)]] <- ccstatus
    phenoOut0[[paste0("X", .phecode)]] <- ccstatus0

    pheinfo[phecode == .phecode, `:=` (
        ncontrols = length(which(ccstatus == 0)),
        ncases = length(which(ccstatus == 1))
    )]
    setTxtProgressBar(pb, i)
} 

fwrite(
    pheinfo,
    glue("data/private/mgi_pheinfo_{version}.txt"),
    sep="\t",row.names=F,col.names=T,quote=T
)

qsave(
    phenoOut,
    glue("data/private/mgi_pim_{version}.qs")
)

qsave(
    phenoOut0,
    glue("data/private/mgi_pim0_{version}.qs")
)

### PHECODE X ------------------------------------------------------------------
phecodesx <- PheWAS::mapCodesToPhecodes(
    input = diag,
    vocabulary.map = ms::phecodex_icdcm_map,
    rollup.map = ms::phecodex_rollup
)

restrictPhecodesxBySex_mod <- function(phenotypes, id.sex, by_var = "id", sex_var = "sex") {
    data <- merge.data.table(
        as.data.table(phenotypes),
        as.data.table(id.sex),
        by = by_var, all.x = TRUE
    )
    # Get the restrictions found in the phenotypes data frame
    current_sex_restriction <- PheWAS::sex_restriction[PheWAS::sex_restriction$phecode %in% unique(data[, phecode]), ] |>
        as.data.table()
    # Get male and female-only phenotypes
    male_only <- current_sex_restriction[current_sex_restriction$male_only, phecode]
    female_only <- current_sex_restriction[current_sex_restriction$female_only, phecode]
    # Set row column matches to NA where inds of a sex meet restricted phenotypes
    data[phecode %in% male_only & sex == "F", phecode := NA]
    data[phecode %in% female_only & sex == "M", phecode := NA]

    na.omit(data)[, (sex_var) := NULL][]
}

phecodesx_restricted <- restrictPhecodesxBySex_mod(
    phenotypes = phecodesx,
    id.sex = demo
)

first_restrictedx <- phecodesx_restricted[phecodesx_restricted[, .I[which.min(dsb)], by = c("id", "phecode")][["V1"]]]

qsave(
    phecodesx_restricted,
    glue("data/private/mgi_phe_dsbx_{version}.qs"),
    nthreads = 4
)

qsave(
    first_restrictedx,
    glue("data/private/mgi_first_phe_dsbx_{version}.qs"),
    nthreads = 4
)

pim0x <- dcast(
    first_restrictedx[, .(id, phecode)],
    id ~ phecode,
    value.var = "phecode",
    fun.aggregate = length,
    fill = 0
)
# collect integer variable names from data.table
int_vars <- names(pim0x)[which(sapply(pim0x, is.integer))]
pim0x[, (int_vars) := lapply(.SD, \(x) as.numeric(x > 0)), .SDcols = int_vars]

qsave(
    pim0x,
    glue("data/private/mgi_pim0x_{version}.qs")
)

message("Done! :)")
