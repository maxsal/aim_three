ms::libri(ms, qs, tidyverse, data.table, cli, glue, survey, mice, patchwork)

sample_sizes <- c(1000, 2500, 5000, 10000)
iterations <- 1000
outcome <- "glucose"
exposure <- "bmi"

res <- map(
    sample_sizes,
    \(i) {
        fread(glue("results/{outcome}_{exposure}_n{i}_i{iterations}_biased_simul_diag.csv"))[, `:=` (size = i)]
    }
)

cov_rate_plot <- rbindlist(res) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "woPRS-imputed",
            method == "Imputed w/ PRS" ~ "PRS-imputed"
        ),
        adj = case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        )
    ) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Covariate-adjusted")),
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))) |>
    filter(adj == "Covariate-adjusted") |>
    dplyr::select(size, adj, method, miss_data, cov_rate, weight) |>
    ggplot(aes(x = size, y = cov_rate, shape = method, color = miss_data, linetype = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        title = stringr::str_wrap("Coverage rate of 95% CI for BMI coefficient for glucose by missing data mechanism and method and weighting under biased sampling", 75),
        x = "Sample size",
        y = "Coverage rate (%)",
        caption = "Dashed line represents 95% coverage rate"
    ) +
    facet_wrap(~weight) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    glue("results/{outcome}_{exposure}_i{iterations}_coverage_rate_biased.pdf"),
    plot = cov_rate_plot,
    width = 8,
    height = 6,
    device = cairo_pdf
)

per_bias_plot <- rbindlist(res) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "woPRS-imputed",
            method == "Imputed w/ PRS" ~ "PRS-imputed"
        ),
        adj = case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        )
    ) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Covariate-adjusted")),
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))) |>
    filter(adj == "Covariate-adjusted") |>
    dplyr::select(size, adj, method, miss_data, per_bias, weight) |>
    ggplot(aes(x = size, y = per_bias/100, shape = method, color = miss_data, linetype = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = stringr::str_wrap("Percent bias for BMI coefficient for glucose by missing data mechanism and method and weighting under biased sampling", 75),
        x = "Sample size",
        y = "Percent bias (%)",
        caption = "Dashed line represents no bias"
    ) +
    facet_wrap(~ weight) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    glue("results/{outcome}_{exposure}_i{iterations}_percent_bias_biased.pdf"),
    plot = per_bias_plot,
    width = 8,
    height = 6,
    device = cairo_pdf
)

leg <- ggpubr::get_legend(cov_rate_plot)

cov_per_comb_plot <- ((cov_rate_plot +
    labs(title = "A. Coverage rate") +
    theme(legend.position = "none"))/
    (per_bias_plot +
    labs(title = "B. Percent bias") +
    theme(legend.position = "none")) /
    leg) +
    # plot_annotation(guides = "collect") +
    plot_layout(heights = c(5, 5, 1)) 

ggsave(
    plot = cov_per_comb_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_coverage_rate_percent_bias_biased.pdf"),
    width = 7,
    height = 6,
    device = cairo_pdf
)

rbindlist(res) |>
    pivot_longer(
        cols = -c(size, adj, method, miss_data, weight),
        names_to = "metric",
        values_to = "value"
    ) |>
    filter(metric != "raw_bias") |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "OD-imputed",
            method == "Imputed w/ PRS" ~ "OD-PRS-imputed"
        )
    ) |>
    dplyr::mutate(
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        metric = factor(metric, levels = c("per_bias", "cov_rate", "avg_wide", "rmse")),
        method = factor(method, levels = c("Complete case", "OD-imputed", "OD-PRS-imputed"))
    ) |>
    arrange(weight, adj, miss_data, metric, method) |>
    pivot_wider(
        names_from = size,
        values_from = value
    ) |>
    dplyr::select(weight, adj, miss_data, metric, method, everything()) |>
    fwrite(glue("results/{outcome}_{exposure}_i{iterations}_summary_biased.csv"))


###########
rbindlist(res) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "woPRS-imputed",
            method == "Imputed w/ PRS" ~ "PRS-imputed"
        ),
        adj = case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        )
    ) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Covariate-adjusted")),
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))) |>
    filter(adj == "Covariate-adjusted") |>
    dplyr::select(size, adj, method, miss_data, avg_wide, weight) |>
    ggplot(aes(x = size, y = avg_wide, shape = method, color = miss_data, linetype = method)) +
    # geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    # scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        title = stringr::str_wrap("Average 95% CI width for BMI coefficient for glucose by missing data mechanism and method and weighting under biased sampling", 75),
        x = "Sample size",
        y = "Average 95% CI width"
    ) +
    facet_wrap(~ weight) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )

ggsave(
    glue("results/{outcome}_{exposure}_i{iterations}_avg_wide_biased.pdf"),
    width = 7,
    height = 5,
    device = cairo_pdf
)

rbindlist(res) |>
    dplyr::mutate(
        method = case_when(
            method == "Complete case" ~ "Complete case",
            method == "Imputed" ~ "woPRS-imputed",
            method == "Imputed w/ PRS" ~ "PRS-imputed"
        ),
        adj = case_when(
            adj == "Unadjusted" ~ "Unadjusted",
            adj == "Adjusted" ~ "Covariate-adjusted"
        )
    ) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Covariate-adjusted")),
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))) |>
    filter(adj == "Covariate-adjusted") |>
    dplyr::select(size, adj, method, miss_data, rmse, weight) |>
    ggplot(aes(x = size, y = rmse, shape = method, color = miss_data, linetype = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    # scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        title = stringr::str_wrap("RMSE for BMI coefficient for glucose by missing data mechanism and method and weighting under biased sampling", 75),
        x = "Sample size",
        y = "RMSE",
        caption = "Dashed line represents no error"
    ) +
    facet_wrap(~ weight) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    glue("results/{outcome}_{exposure}_i{iterations}_rmse_biased.pdf"),
    width = 7,
    height = 5,
    device = cairo_pdf
)
