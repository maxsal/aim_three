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

ex_prs_files <- fread("data/public/ex_prs_files.csv")

these_ex_prs <- ex_prs_files |>
    filter(use == 1) |>
    pull(trait)

ex_prs <- map(
    1:nrow(ex_prs_files),
    \(r) {
        tmp <- fread(ex_prs_files$prs_file[r])
        names(tmp) <- c("id", janitor::make_clean_names(ex_prs_files$trait[r]))
        # names(tmp) <- c("id", ex_prs_files$label[r])
        return(tmp)
    }
) |>
    reduce(full_join, by = "id") |>
    select(any_of(c("id", these_ex_prs)))

lung_can_prs <- fread(paste0(prs_path, "PRSWEB_PHECODE165.1_LUNG-CANCER-MESOT_PRS-CS_MGI_20200608_PRS.txt"))

ex_prs_cor <- cor(ex_prs[, !c("id")])

library(ggcorrplot)

ggcorrplot(
    ex_prs_cor,
    hc.order = TRUE,
    # outline.col = "white",
    show.diag = FALSE
) +
    labs(
        title = "ExPRS v. ExPRS"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank()
    )
