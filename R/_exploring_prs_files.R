ms::libri(ms, qs, data.table, tidyverse, glue, cli)

# paths
base_path <-  "/net/junglebook/home/larsf/students/shareMax"
prs_path  <- paste0(base_path, "/PRS/")
exprs_path <- paste0(base_path, "/ExPRS/")

# exprs
traits <- c("Vitamin-D", "Vitamin-B12", "Triglycerides", "TotalCholesterol",
    "Thinness", "T2d", "Systolic-blood-pressure", "Smoker-Yes", "LDL",
    "insomnia", "Hypertension", "Height", "HDL", "Glucose",
    "FPG_MAGIC-FastingGlucose", "FPG_FGovertime-corrected",
    "Estradioal-Oestradiol", "eGFR", "Alcohol-Amount")

ex_trait <- "Vitamin-D"
prefix     <- "ExPRSweb_"
approaches <- c("PT", "PRSCS", "PLINK", "LASSOSUM", "DBSLMM")
ex_approach <- "PRSCS"
freeze     <- "3"
date       <- "20211120"

other <- "100021-raw"

filename <- paste0(prefix, ex_trait, "_", other, "_", ex_approach, "_MGI-Freeze", freeze, "_", date, "_PRS.txt")

tmp <- fread(paste0(exprs_path, filename))
