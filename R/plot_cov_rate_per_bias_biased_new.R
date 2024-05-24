# libraries --------------------------------------------------------------------
ms::libri(ms, qs, tidyverse, data.table, cli, glue, survey, mice, patchwork)

# settings
sample_sizes <- c(1000, 2500, 5000, 10000)
iterations <- 1000
outcome <- "glucose"
exposure <- "bmi"

# data -------------------------------------------------------------------------
## exposure missing only
ex_only_res <- map(
    sample_sizes,
    \(i) {
        fread(glue("results/simulations/exposure_only/{outcome}_{exposure}_n{i}_i{iterations}_biased_simul_diag.csv"))[, `:=` (size = i)]
    }
)

## exposure and outcome missing 
ex_out_res <- map(
    sample_sizes,
    \(i) {
        fread(glue("results/simulations/exposure_and_outcome/{outcome}_{exposure}_n{i}_i{iterations}_exoutmiss_biased_simul_diag.csv"))[, `:=`(size = i)]
    }
)

# plots ------------------------------------------------------------------------
## coverage rate plot
make_diag_plot <- function(
    results_list,
    diagnostic_var,
    title = NULL,
    x_lab = NULL,
    y_lab = NULL,
    caption = NULL,
    title_wrap = 75,
    hline = NULL
) {
    rbindlist(results_list) |>
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
            miss_dat = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
            method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))
        ) |>
        filter(adj == "Covariate-adjusted") |>
        dplyr::select(size, adj, method, miss_dat, weight, {{ diagnostic_var }}) |>
        ggplot(aes(x = size, y = {{ diagnostic_var }}, shape = method, color = miss_dat, linetype = method)) +
        geom_hline(yintercept = hline, linetype = "dashed") +
        geom_line(linewidth = 1) +
        geom_point(size = 4) +
        scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        facet_wrap(~weight) +
        labs(
            title = stringr::str_wrap(title, title_wrap),
            x = x_lab,
            y = y_lab,
            caption = caption
        ) +
        scale_color_ms() +
        theme_ms() +
        theme(
            legend.title = element_blank()
        )
}


### exposure only
ex_only_cov_rate_plot <- make_diag_plot(
    ex_only_res,
    cov_rate,
    title = "Coverage rate for BMI coefficient for glucose by missing data mechanism and method under random sampling",
    x_lab = "Sample size",
    y_lab = "Coverage rate (%)",
    caption = "Dashed line represents nominal coverage rate",
    hline = 0.95
)
ggsave(
    plot = ex_only_cov_rate_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_exmiss_coverage_rate_biased.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)


### exposure and outcome
ex_and_out_cov_rate_plot <- make_diag_plot(
    ex_out_res,
    cov_rate,
    title = "Coverage rate for BMI coefficient for glucose by missing data mechanism and method under random sampling",
    x_lab = "Sample size",
    y_lab = "Coverage rate (%)",
    caption = "Dashed line represents nominal coverage rate",
    hline = 0.95
)
ggsave(
    plot = ex_and_out_cov_rate_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_exoutmiss_coverage_rate_biased.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)


# percent bias plot
make_per_bias_plot <- function(
    results_list
) {
    rbindlist(results_list) |>
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
            miss_dat = factor(miss_data, levels = c("MCAR", "MAR", "MNAR")),
            method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))
        ) |>
        filter(adj == "Covariate-adjusted") |>
        dplyr::select(size, adj, method, miss_dat, weight, per_bias) |>
        ggplot(aes(x = size, y = per_bias / 100, shape = method, color = miss_dat, linetype = method)) +
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
        facet_wrap(~ weight) +
        scale_color_ms() +
        theme_ms() +
        theme(
            legend.title = element_blank()
        )
}

## exposure only
ex_only_per_bias_plot <- make_per_bias_plot(ex_only_res)
ggsave(
    plot = ex_only_per_bias_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_exmiss_percent_bias_biased.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)

## exposure and outcome
ex_out_per_bias_plot <- make_per_bias_plot(ex_out_res)
ggsave(
    plot = ex_out_per_bias_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_exoutmiss_percent_bias_biased.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)

# coverage rate and percent bias combined plot

## exposure only
ex_leg <- ggpubr::get_legend(ex_only_cov_rate_plot)

ex_only_cov_per_comb_plot <- ((ex_only_cov_rate_plot +
    labs(title = "A. Coverage rate") +
    theme(legend.position = "none")) /
    (ex_only_per_bias_plot +
    labs(title = "B. Percent bias") +
    theme(legend.position = "none")) /
    ex_leg) +
    plot_layout(heights = c(5, 5, 1))
ggsave(
    plot = ex_only_cov_per_comb_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_exmiss_coverage_rate_percent_bias_biased.pdf"),
    width = 7,
    height = 6,
    device = cairo_pdf
)


## exposure and outcome
ex_out_leg <- ggpubr::get_legend(ex_and_out_cov_rate_plot)

ex_out_cov_per_comb_plot <- ((ex_and_out_cov_rate_plot +
    labs(title = "A. Coverage rate") +
    theme(legend.position = "none")) /
    (ex_out_per_bias_plot +
    labs(title = "B. Percent bias") +
    theme(legend.position = "none")) /
    ex_out_leg) +
    plot_layout(heights = c(5, 5, 1))
ggsave(
    plot = ex_out_cov_per_comb_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_exoutmiss_coverage_rate_percent_bias_biased.pdf"),
    width = 7,
    height = 6,
    device = cairo_pdf
)

layout <- "
AAABBB
CCCDDD
##EE##
"

cov_per_comb_plot <- (ex_only_cov_rate_plot +
    labs(title = "A. Exposure only, coverage rate", caption = "") +
    theme(legend.position = "none")) + (ex_and_out_cov_rate_plot +
    labs(title = "B. Exposure and outcome, coverage rate", caption = "") +
    theme(legend.position = "none")) +

    (ex_only_per_bias_plot +
        labs(title = "C. Exposure only, percent bias", caption = "") +
        scale_y_continuous(limits = c(0, 1.6), labels = scales::percent) +
        # ylim(0, 0.8) +
        theme(legend.position = "none")) + (ex_out_per_bias_plot +
        scale_y_continuous(limits = c(0, 1.6), labels = scales::percent) +
        labs(title = "D. Exposure and outcome, percent bias", caption = "") +
        theme(legend.position = "none")) +
    ex_out_leg +
    plot_layout(design = layout, heights = c(10, 10, 1))

ggsave(
    plot = cov_per_comb_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_combined_coverage_rate_percent_bias_biased.pdf"),
    width = 12,
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

#########
avg_wide_plot <- rbindlist(res) |>
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
        miss_dat = factor(miss_dat, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))) |>
    filter(adj == "Covariate-adjusted") |>
    dplyr::select(size, adj, method, miss_dat, avg_wide) |>
    ggplot(aes(x = size, y = avg_wide, shape = method, color = miss_dat, linetype = method)) +
    # geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    # scale_y_continuous(labels = scales::percent) +
    labs(
        title = stringr::str_wrap("Average 95% CI width for BMI coefficient for glucose by missing data mechanism and method under random sampling", 75),
        x = "Sample size",
        y = "Average 95% CI width"
    ) +
    # facet_wrap(~adj) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )

ggsave(
    plot = avg_wide_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_avg_wide_random.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)

rmse_plot <- rbindlist(res) |>
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
        miss_dat = factor(miss_dat, levels = c("MCAR", "MAR", "MNAR")),
        method = factor(method, levels = c("Complete case", "woPRS-imputed", "PRS-imputed"))) |>
    filter(adj == "Covariate-adjusted") |>
    dplyr::select(size, adj, method, miss_dat, rmse) |>
    ggplot(aes(x = size, y = rmse, shape = method, color = miss_dat, linetype = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = sample_sizes, limits = c(0, max(sample_sizes)), labels = scales::comma) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
        title = stringr::str_wrap("RMSE for BMI coefficient for glucose by missing data mechanism and method under random sampling", 75),
        x = "Sample size",
        y = "RMSE",
        caption = "Dashed line represents no error"
    ) +
    # facet_wrap(~adj) +
    scale_color_ms() +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )

ggsave(
    plot = rmse_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_rmse_random.pdf"),
    width = 8,
    height = 5,
    device = cairo_pdf
)

leg <- ggpubr::get_legend(avg_wide_plot)

wide_rmse_comb_plot <- ((avg_wide_plot +
    labs(title = "A. Average 95% CI width") +
    theme(legend.position = "none")) /
    (rmse_plot +
        labs(title = "B. RMSE") +
        theme(legend.position = "none")) /
    leg) +
    plot_layout(heights = c(5, 5, 1))

ggsave(
    plot = wide_rmse_comb_plot,
    glue("results/{outcome}_{exposure}_i{iterations}_avg_wide_rmse_random.pdf"),
    width = 8,
    height = 7,
    device = cairo_pdf
)
