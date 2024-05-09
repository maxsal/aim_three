ms::libri(ms, qs, tidyverse, data.table, cli, glue, survey, mice, patchwork)

sample_sizes <- c(1000, 2500, 5000, 10000)
iterations <- 1000
outcome <- "glucose"
exposure <- "bmi"

res <- map(
    sample_sizes,
    \(i) {
        fread(glue("results/{outcome}_{exposure}_n{i}_i{iterations}_diag_random.csv"))[, `:=` (size = i)]
    }
)

cov_rate_plot <- rbindlist(res) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "OD-imputed",
            method == "Imputed w/ PRS" ~ "OD-PRS-imputed"
        ),
        adj = case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        )
    ) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Covariate-adjusted")),
        miss_dat = factor(miss_dat, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "OD-imputed", "OD-PRS-imputed"))) |>
    dplyr::select(size, adj, method, miss_dat, cov_rate) |>
    ggplot(aes(x = size, y = cov_rate, shape = method, color = miss_dat, linetype = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        title = stringr::str_wrap("Coverage rate of 95% CI for BMI coefficient for glucose by missing data mechanism and method under random sampling", 75),
        x = "Sample size",
        y = "Coverage rate (%)",
        caption = "Dashed line represents 95% coverage rate"
    ) +
    facet_wrap(~adj) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    plot = cov_rate_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_coverage_rate_random.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)

per_bias_plot <- rbindlist(res) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "OD-imputed",
            method == "Imputed w/ PRS" ~ "OD-PRS-imputed"
        ),
        adj = case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        )
    ) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Covariate-adjusted")),
        miss_dat = factor(miss_dat, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "OD-imputed", "OD-PRS-imputed"))
    ) |>
    dplyr::select(size, adj, method, miss_dat, per_bias) |>
    ggplot(aes(x = size, y = per_bias/100, shape = method, color = miss_dat, linetype = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = stringr::str_wrap("Percent bias for BMI coefficient for glucose by missing data mechanism and method under random sampling", 75),
        x = "Sample size",
        y = "Percent bias (%)",
        caption = "Dashed line represents no bias"
    ) +
    facet_wrap(~ adj) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    plot = per_bias_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_percent_bias_random.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)

##

leg <- cowplot::get_legend(cov_rate_plot)

cov_per_comb_plot <- ((cov_rate_plot +
    labs(title = "A. Coverage rate") +
    theme(legend.position = "none"))/
    (per_bias_plot +
    labs(title = "B. Percent bias") +
    theme(legend.position = "none")) /
    leg) +
    plot_layout(heights = c(5, 5, 1))

ggsave(
    plot = cov_per_comb_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_coverage_rate_percent_bias_random.pdf"),
    width = 8,
    height = 7,
    device = cairo_pdf
)

##

rbindlist(res) |>
    select(adj, miss_dat, method, per_bias, size) |>
    pivot_wider(
        names_from = size,
        values_from = per_bias
    )

rbindlist(res) |>
    pivot_longer(
        cols = -c(size, adj, method, miss_dat),
        names_to = "metric",
        values_to = "value"
    ) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "OD-imputed",
            method == "Imputed w/ PRS" ~ "OD-PRS-imputed"
        )
    ) |>
    filter(metric != "raw_bias") |>
    dplyr::mutate(
        miss_dat = factor(miss_dat, levels = c("MCAR", "MAR", "MNAR")),
        metric = factor(metric, levels = c("per_bias", "cov_rate", "avg_wide", "rmse")),
        method = factor(method, levels = c("Complete case", "OD-mputed", "OD-PRS-imputed"))
    ) |>
    arrange(adj, miss_dat, metric, method) |>
    pivot_wider(
        names_from = size,
        values_from = value
    ) |>
    dplyr::select(adj, miss_dat, metric, method, everything()) |>
    fwrite(glue("results/{outcome}_{exposure}_i{iterations}_summary_random.csv"))
