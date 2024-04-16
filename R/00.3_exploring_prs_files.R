# libraries --------------------------------------------------------------------
ms::libri(ms, qs, data.table, tidyverse, glue, cli, ggcorrplot)

# paths
base_path  <- "/net/junglebook/home/larsf/students/shareMax"
prs_path   <- paste0(base_path, "/PRS/")
exprs_path <- paste0(base_path, "/ExPRS/")

# data -------------------------------------------------------------------------
ex_prs_files <- fread("data/public/ex_prs_files.csv")
prsweb_files <- fread("data/public/prsweb_files.csv", colClasses = "character") |>
    mutate(
        file = paste0(prs_path, prsweb_file)
    ) |>
    filter(prsweb_file != "")

these_ex_prs <- ex_prs_files |>
    filter(use == 1) |>
    pull(trait)
these_prsweb <- prsweb_files |>
    filter(prsweb_file != "") |>
    pull(cancer) |>
    janitor::make_clean_names()

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

qsave(
    ex_prs,
    file = "data/private/ex_prs_processed.qs"
)

################################################################################
### CANCER PRSWEB FILES HAVE OLD IDS - NO WAY OF LINKING TO EXPRS FILES ########
################################################################################
# pull cancer prsweb files
# prsweb <- map(
#     1:nrow(prsweb_files),
#     \(r) {
#         tmp <- fread(prsweb_files$file[r])
#         names(tmp) <- c("id", janitor::make_clean_names(prsweb_files$cancer[r]))
#         return(tmp)
#     }
# ) |>
#     reduce(full_join, by = "id") |>
#     select(any_of(c("id", these_prsweb))
# )

# # combine prsweb and ex_prs ----------------------------------------------------
# prs <- ex_prs |>
#     full_join(prsweb, by = "id")

################################################################################
# # correlation matrix ---------------------------------------------------------
################################################################################
# quant_traits_from_paper <- rev(c("alc_amount", "sbp_raw", "dbp_raw", "height",
#     "bmi", "glucose", "vit_b12", "vit_d_raw", "creatinine_qntl", "egfr",
#     "crp_qntl", "chol_qntl", "tri_qntl", "ldl", "hdl"))

# var_labels <- ex_prs_files |>
#     filter(trait %in% quant_traits_from_paper) |>
#     select(label, trait) |>
#     deframe()
# ex_prs_cor <- cor(ex_prs |>
#     select(any_of(quant_traits_from_paper)) |>
#     rename(!!!var_labels)
# )

# ggcorrplot(
#     ex_prs_cor,
#     colors = c("#0072B2", "white", "#D55E00"),
#     show.diag = FALSE
# ) +
#     labs(
#         title = "ExPRS v. ExPRS"
#     ) +
#     theme(
#         plot.title = element_text(hjust = 0.5),
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold"),
#         panel.grid.major = element_blank()
#     )
