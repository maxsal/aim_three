ms::libri(ms, qs, tidyverse, data.table, cli, glue, survey, mice, patchwork)

sample_sizes <- c(1000, 2500, 5000, 10000)
iterations <- 1000
outcome <- "glucose"
exposure <- "bmi"

res <- map(
    sample_sizes,
    \(i) {
        fread(glue("data/public/{outcome}_{exposure}_n{i}_i{iterations}_diag.csv"))[, `:=` (size = i)]
    }
)

rbindlist(res) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Adjusted")),
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "Imputed", "Imputed w/ PRS"))) |>
    dplyr::select(size, adj, method, miss_data, cov_rate, weight) |>
    ggplot(aes(x = size, y = cov_rate, shape = method, color = miss_data, linetype = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        title = stringr::str_wrap("Coverage rate of 95% CI for BMI coefficient for glucose by missing data mechanism and method and weighting under biased sampling", 75),
        x = "Sample size",
        y = "Coverage rate (%)",
        caption = "Dashed line represents 95% coverage rate"
    ) +
    facet_wrap(~adj+weight) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    glue("results/{outcome}_{exposure}_i{iterations}_coverage_rate_biased.pdf"),
    width = 8,
    height = 6,
    device = cairo_pdf
)

rbindlist(res) |>
    dplyr::mutate(
        adj = factor(adj, levels = c("Unadjusted", "Adjusted")),
        miss_data = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "Imputed", "Imputed w/ PRS"))
    ) |>
    dplyr::select(size, adj, method, miss_data, per_bias, weight) |>
    ggplot(aes(x = size, y = per_bias/100, shape = method, color = miss_data, linetype = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = stringr::str_wrap("Percent bias for BMI coefficient for glucose by missing data mechanism and method and weighting under biased sampling", 75),
        x = "Sample size",
        y = "Percent bias (%)",
        caption = "Dashed line represents no bias"
    ) +
    facet_wrap(~ adj + weight) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )
ggsave(
    glue("results/{outcome}_{exposure}_i{iterations}_percent_bias_biased.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)
