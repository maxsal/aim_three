ms::libri(ms, qs, tidyverse, data.table, cli, glue, survey, mice, patchwork)

sample_sizes <- c(1000, 2500, 5000, 10000)
iterations <- 1000
outcome <- "glucose"
exposure <- "bmi"



res <- map(
    sample_sizes,
    \(i) {
        fread(glue("results/{outcome}_{exposure}_n{i}_i{iterations}_biased_simul_res.csv"))[, `:=` (size = i)]
    }
) |> set_names(sample_sizes)


fix_for_plot <- function(x) {
    full_un  <- x[method == "Full" & adj == "Unadjusted", ]
    full_adj <- x[method == "Full" & adj == "Adjusted", ]

    full_un2 <- bind_rows(
        full_un |> dplyr::mutate(weight = "Unweighted", miss_data = "Full"),
        full_un |> dplyr::mutate(weight = "Weighted", miss_data = "Full")
    )
    full_adj2 <- bind_rows(
        full_adj |> dplyr::mutate(weight = "Unweighted", miss_data = "Full"),
        full_adj |> dplyr::mutate(weight = "Weighted", miss_data = "Full")
    )

    out <- x[method != "Full", ]
    bind_rows(
        full_un2,
        full_adj2,
        out
    )
}

plots <- map(
    names(res),
    \(s) {
        fix_for_plot(res[[s]]) |>
        mutate(
            method = factor(method, c("Full", "Complete case", "Imputed", "Imputed w/ PRS")),
            miss_data = factor(miss_data, c("Full", "MCAR", "MAR", "MNAR"))
        ) |>
        ggplot(aes(x = est, y = ..scaled.., color = miss_data, linetype = method)) +
        geom_density(linewidth = 1) +
        facet_wrap(~ adj + weight) +
        scale_color_ms(black_first = TRUE) +
        labs(
            title = glue("{s} samples"),
            x = "Estimate",
            y = "Density (scaled)"
        ) +
        theme_ms() +
        theme(
            legend.title = element_blank()
        )
    }

)

biased_simul_plots <- plots |> patchwork::wrap_plots() +
    plot_layout(guides = "collect") &
    theme(
        legend.position = "bottom"
    )
ggsave(
    "results/biased_simul_plots.pdf",
    biased_simul_plots,
    width = 12,
    height = 8,
    device = cairo_pdf
)
