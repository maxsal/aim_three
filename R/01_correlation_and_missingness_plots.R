# libraries --------------------------------------------------------------------
ms::libri(ms, qs, data.table, tidyverse, ggcorrplot, visdat, naniar,
    labelled, patchwork, mice)

# data -------------------------------------------------------------------------
demo         <- qread("data/private/mgi/20220822/demo_com_labs_ehr_20220822.qs")
ex_prs       <- qread("data/private/ex_prs_processed.qs")
ex_prs_files <- fread("data/public/ex_prs_files.csv")

both_ids <- intersect(demo[, id], ex_prs[, id])

demo   <- demo |> filter(id %in% both_ids)
ex_prs <- ex_prs |> filter(id %in% both_ids)

demo_var_labels <- fread("data/public/demo_var_labels.csv")

these_demo_vars <- demo_var_labels |>
    dplyr::filter(quant == 1) |>
    dplyr::pull(demo_var) |>
    rev()
these_exprs_vars <- rev(c(
    "alc_amount", "sbp_raw", "dbp_raw", "height",
    "bmi", "glucose", "vit_b12", "vit_d_raw", "creatinine_qntl", "egfr",
    "crp_qntl", "chol_qntl", "tri_qntl", "ldl", "hdl"
))

demo_var_labs <- demo_var_labels |>
    select(label, demo_var) |>
    deframe()
these_demo_var_labs <- demo_var_labels |>
    filter(quant == 1) |>
    select(label, demo_var) |>
    deframe()
exprs_var_labs <- ex_prs_files |>
    filter(trait %in% these_exprs_vars) |>
    select(label, trait) |>
    deframe()

# correlation plots ------------------------------------------------------------
demo_cor_plot <- demo |>
    dplyr::select(all_of(these_demo_vars)) |>
    dplyr::rename(!!!these_demo_var_labs) |>
    cor(use = "pairwise.complete.obs") |>
    ggcorrplot(
        colors = c("#0072B2", "white", "#D55E00"),
        show.diag = FALSE
    ) +
    labs(
        title = "Exposure v. Exposure"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank()
    )

ex_prs_cor_plot <- ex_prs |>
    dplyr::select(all_of(these_exprs_vars)) |>
    dplyr::rename(!!!exprs_var_labs) |>
    cor(use = "pairwise.complete.obs") |>
    ggcorrplot(
        colors = c("#0072B2", "white", "#D55E00"),
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

ex_prs_demo_cor_plot <- cor(
    ex_prs |>
        dplyr::select(all_of(these_exprs_vars)) |>
        dplyr::rename(!!!exprs_var_labs),
    demo |>
        dplyr::select(all_of(these_demo_vars)) |>
        dplyr::rename(!!!these_demo_var_labs),
    use = "pairwise.complete.obs"
) |>
    ggcorrplot(
        colors = c("#0072B2", "white", "#D55E00")
    ) +
    labs(
        title = "ExPRS v. Exposure"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank()
    )

patched <- (
    (demo_cor_plot +
        labs(title = "A. Exposure v. Exposure") +
        theme(
            plot.title = element_text(face = "bold", hjust = 0),
            legend.position = "none"
        )) +
    (ex_prs_cor_plot +
        labs(title = "B. ExPRS v. ExPRS") +
        theme(plot.title = element_text(face = "bold", hjust = 0))) +
    (ex_prs_demo_cor_plot +
        labs(title = "C. ExPRS v. Exposure") +
        theme(plot.title = element_text(face = "bold", hjust = 0),
            legend.position = "none"
        ))
    ) +
    plot_layout(nrow = 1)

ggsave(
    "results/correlation_plots.pdf",
    patched,
    width = 14,
    height = 5.5
)

# missingness plot -------------------------------------------------------------
cols <- c(
    "None" = "#009E73",
    "Some" = "black"
)
missing_plot <- demo |>
    rename(!!!demo_var_labs) |>
    map_dbl(\(x) sum(is.na(x)) / length(x)) |>
    enframe() |>
    mutate(not_missing = ifelse(value == 0, "None", "Some")) |>
    ggplot(aes(x = reorder(name, value), y = value)) +
    geom_segment(aes(xend = name, y = 0, yend = value)) +
    geom_point(aes(color = not_missing)) +
    coord_flip() +
    scale_color_manual(values = cols) +
    labs(
        title = "Proportion missing by variable",
        # subtitle = paste0("n = ", format(nrow(demo), big.mark = ",")),
        x = "Variable",
        y = "Proportion Missing"
    ) +
    ms::theme_ms() +
    theme(
        panel.grid.major.y = element_blank(),
        legend.position = "none"
    )

ggsave(
    "results/missingness_plot.pdf",
    missing_plot,
    width = 6,
    height = 7.5
)

# save analytic data -----------------------------------------------------------
old_ex_prs_names        <- names(ex_prs)[names(ex_prs) != "id"]
names(old_ex_prs_names) <- paste0(old_ex_prs_names, "_prs")

analytic_full <- full_join(
    demo,
    ex_prs |>
        rename(!!!old_ex_prs_names),
    by = "id"
)

keep_these_vars <- demo |>
    map_dbl(\(x) sum(is.na(x)) / length(x)) |>
    enframe() |>
    filter(value <= 0.25) |>
    pull(name)

analytic_sub <- analytic_full |>
    select(any_of(keep_these_vars)) |>
    na.omit()

qsave(
    analytic_full,
    "data/private/analytic_full.qs"
)

qsave(
    analytic_sub,
    "data/private/analytic_sub.qs"
)
