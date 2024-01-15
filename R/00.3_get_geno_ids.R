# Script for identifying set of MGI samples with ExPRS data
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

exprs_file <- "/net/junglebook/home/larsf/students/shareMax/ExPRS/ExPRSweb_T2D_29358691-GCST005413_DBSLMM_MGI-Freeze3_20211120_PRS.txt"
freeze4_file <- "/net/junglebook/magic_data/MGI_GenotypeData_Freeze4/MGI_Freeze4_Sample_Info_60215samples.txt"

exprs <- fread(exprs_file)
setnames(exprs, "IID", "id")
freeze4 <- fread(freeze4_file)
setnames(freeze4, "Deid_ID", "id")

geno_ids <- data.table(
    id = intersect(
    exprs[, id],
    freeze4[, id]
    )
)

fwrite(
    geno_ids,
    glue("data/private/mgi_geno_ids_{version}.txt")
)

message("Done! :)")
