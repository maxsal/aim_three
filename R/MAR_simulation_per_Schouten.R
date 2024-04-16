library(mice)
library(tidyverse)
library(MASS)
library(furrr)
plan(multisession, workers = 6)

rho <- 0.2

sample <- MASS::mvrnorm(
  n = 1000,
  mu = c(5, 5, 10),
  Sigma = matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), ncol = 3)
) |> as.data.frame()
names(sample) <- c("Y1", "Y2", "X1")



pattern <- c(0, 0, 1)
weights <- c(0, 0, 1)

mars <- map(
  1:100,
  \(i) {
    ampute(
      data = sample, 
      prop = 0.3,
      patterns = pattern,
      mech = "MAR",
      weights = weights
    )
  },
  .progress = TRUE
)

imps <- future_map(
  mars,
  \(i) {
    mice(i$amp, printFlag = FALSE)
  },
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

# mean
true_mean <- mean(sample[["Y1"]])

cc_mean_est <- map_dbl(
  mars,
  \(i) {
    mean(i$amp[, 1], na.rm = TRUE)
  }
)
imp_mean_est <- map_dbl(
  imps,
  \(x) {
    mean(map_dbl(seq_len(x$m), \(i) mean(complete(x, i)[, 1])))
  }
)

# Y1 ~ X1 association
true_un_ass  <- glm(Y1 ~ X1, data = sample)$coefficients[["X1"]]
true_adj_ass <- glm(Y1 ~ X1 + Y2, data = sample)$coefficients[["X1"]]

cc_un_ass <- map_dbl(
  mars,
  \(i) {
    glm(Y1 ~ X1, data = i$amp)$coefficients[["X1"]]
  }
)
cc_adj_ass <- map_dbl(
  mars,
  \(i) {
    glm(Y1 ~ X1 + Y2, data = i$amp)$coefficients[["X1"]]
  }
)

imp_un_ass <- map_dbl(
  imps,
  \(x) {
    mean(
      map_dbl(
        seq_len(x$m),
        \(i) glm(Y1 ~ X1, data = complete(x, i))$coefficients[["X1"]]
      )
    )
  }
)

imp_adj_ass <- map_dbl(
  imps,
  \(x) {
    mean(
      map_dbl(
        seq_len(x$m),
        \(i) glm(Y1 ~ X1 + Y2, data = complete(x, i))$coefficients[["X1"]]
      )
    )
  }
)

mean_data <- bind_rows(
  tibble(bias_est = cc_mean_est - true_mean, type = "Complete case"),
  tibble(bias_est = imp_mean_est - true_mean, type = "Imputed")
) 
mean_data_mean <- mean_data |>
  summarize(bias_est = mean(bias_est, na.rm = TRUE), .by = c("type"))

mean_data |>
  ggplot(aes(x = bias_est, color = type, linetype = type)) +
  geom_vline(xintercept = 0, linewidth = 1) +
  geom_vline(data = mean_data_mean, aes(xintercept = bias_est, color = type), linewidth = 1, linetype = 2) +
  geom_density(linewidth = 1) +
  labs(
    title = "Bias in Y1 mean estimation under missing scenarios",
    subtitle = c("n = 1000; m = 5; miss_prop = 0.3; amputations = 100"),
    x = "Bias (mean_est - mean_true)",
    y = "Density",
    captione = "Dashed vertical lines represents mean"
  ) +
  ms::scale_color_ms() +
  ms::theme_ms() +
  theme(legend.title = element_blank())
ggsave(
  "~/Downloads/mean_bias_simul.pdf",
  width = 7, height = 5, device = cairo_pdf
)

beta_data <- bind_rows(
  tibble(bias_est = cc_un_ass - true_un_ass, type = "Complete case", cov = "Unadjusted"),
  tibble(bias_est = cc_adj_ass - true_adj_ass, type = "Complete case", cov = "Adjusted"),
  tibble(bias_est = imp_un_ass - true_un_ass, type = "Imputed", cov = "Unadjusted"),
  tibble(bias_est = imp_adj_ass - true_adj_ass, type = "Imputed", cov = "Adjusted")
) 

beta_data_mean <- beta_data |>
  summarize(bias_est = mean(bias_est, na.rm = TRUE), .by = c("type", "cov"))

beta_data |>
  ggplot(aes(x = bias_est, color = type, linetype = type)) +
  geom_vline(xintercept = 0) +
  geom_vline(data = beta_data_mean, aes(xintercept = bias_est, color = type), linewidth = 1, linetype = 2) +
  geom_density(linewidth = 1) +
  facet_wrap(~factor(cov, c("Unadjusted", "Adjusted")), ncol = 1) +
  labs(
    title = "Bias in X1 beta coefficient for Y1 estimation under missing scenarios",
    subtitle = c("n = 1000; m = 5; miss_prop = 0.3; amputations = 100"),
    x = "Bias (Beta_est - Beta_true)",
    y = "Density"
  ) +
  ms::scale_color_ms() +
  ms::theme_ms() +
  theme(legend.title = element_blank())
ggsave(
  "~/Downloads/beta_bias_simul.pdf",
  width = 7, height = 5, device = cairo_pdf
)

## manual missingness
scores <- -1 + scale(sample[["X1"]])[, 1]
miss_probs <- plogis(scores)
miss_ind <- rbinom(length(miss_probs), 1, miss_probs)

man <- sample

mans <- map(
  1:1000,
  \(i) {
    miss_probs <- plogis(scores)
    miss_ind <- rbinom(length(miss_probs), 1, miss_probs)
    man[miss_ind == 1, "Y1"] <- NA
    man[miss_ind == 1, "Y2"] <- NA
    man
  }
)

mnar_scores <- -1 + scale(sample[["X1"]])[, 1] + scale(sample[["Y1"]])[, 1]
mnar_miss_probs <- plogis(mnar_scores)
mnar_miss_ind <- rbinom(length(mnar_miss_probs), 1, mnar_miss_probs)

mnar_man <- sample

mnar_mans <- map(
  1:100,
  \(i) {
    mnar_miss_probs <- plogis(mnar_scores)
    mnar_miss_ind <- rbinom(length(mnar_miss_probs), 1, mnar_miss_probs)
    mnar_man[mnar_miss_ind == 1, "Y1"] <- NA
    mnar_man[mnar_miss_ind == 1, "Y2"] <- NA
    mnar_man
  }
)

cc_man_mean_est <- map_dbl(
  mans,
  \(i) {
    mean(i[["Y1"]], na.rm = TRUE)
  }
)
cc_mnar_man_mean_est <- map_dbl(
  mnar_mans,
  \(i) {
    mean(i[["Y1"]], na.rm = TRUE)
  }
)

man_imp <- map(
  mans,
  \(i) {
    mice(i, printFlag = FALSE)
  },
  .progress = TRUE
)

imp_man_mean_est <- map_dbl(
  man_imp,
  \(x) {
    map_dbl(
      seq_len(x$m),
      \(i) {
        mean(complete(x, i)[["Y1"]])
      }
    ) |> mean()
  }
)


cc_man_un_beta_est <- map_dbl(
  mans,
  \(i) {
    glm(Y1 ~ X1, data = i)$coefficients[["X1"]]
  }
)
cc_man_adj_beta_est <- map_dbl(
  mans,
  \(i) {
    glm(Y1 ~ X1 + Y2, data = i)$coefficients[["X1"]]
  }
)

imp_man_un_beta_est <- map(
  man_imp,
  \(x) {
    fit <- with(x,glm(Y1 ~ X1))
    tab <- summary(pool(fit), "all", conf.int = TRUE) |> as_tibble()
    tab |>
      filter(term == "X1") |>
      dplyr::select(term, est = estimate, lo = `2.5 %`, hi = `97.5 %`)
  }
)


true_un_ass
tmp2 <- bind_rows(imp_man_un_beta_est)

RB <- mean(tmp2$est) - true_un_ass
PB <- 100 * abs((mean(tmp2$est) - true) / true)
CR <- mean(tmp2$lo < true_un_ass & true_un_ass < tmp2$hi)
AW <- mean(tmp2$hi - tmp2$lo)
RMSE <- sqrt(mean(tmp2$est - true_un_ass)^2)

tibble(RB, PB, CR, AW, RMSE)


imp_man_adj_beta_est <- map_dbl(
  man_imp,
  \(x) {
    map_dbl(
      seq_len(x$m),
      \(i) {
        glm(Y1 ~ X1 + Y2, data = complete(x, i))$coefficients[["X1"]]
      }
    ) |> mean()
  }
)

man_mean_data <- bind_rows(
  tibble(bias_est = cc_man_mean_est - true_mean, type = "Complete case"),
  tibble(bias_est = imp_man_mean_est - true_mean, type = "Imputed")
) 
man_mean_data_mean <- man_mean_data |>
  summarize(bias_est = mean(bias_est, na.rm = TRUE), .by = c("type"))

man_mean_data |>
  ggplot(aes(x = bias_est, color = type, linetype = type)) +
  geom_vline(xintercept = 0, linewidth = 1) +
  geom_vline(data = man_mean_data_mean, aes(xintercept = bias_est, color = type), linewidth = 1, linetype = 2) +
  geom_density(linewidth = 1) +
  labs(
    title = "Bias in Y1 mean estimation under missing scenarios",
    subtitle = c("n = 1000; m = 5; miss_prop = 0.3; MAR dataset = 100; manual missingness"),
    x = "Bias (mean_est - mean_true)",
    y = "Density",
    captione = "Dashed vertical lines represents mean"
  ) +
  ms::scale_color_ms() +
  ms::theme_ms() +
  theme(legend.title = element_blank())
ggsave(
  "~/Downloads/man_mean_bias_simul.pdf",
  width = 7, height = 5, device = cairo_pdf
)

man_beta_data <- bind_rows(
  tibble(bias_est = cc_man_un_beta_est - true_un_ass, type = "Complete case", cov = "Unadjusted"),
  tibble(bias_est = cc_man_adj_beta_est - true_adj_ass, type = "Complete case", cov = "Adjusted"),
  tibble(bias_est = imp_man_un_beta_est - true_un_ass, type = "Imputed", cov = "Unadjusted"),
  tibble(bias_est = imp_man_adj_beta_est - true_adj_ass, type = "Imputed", cov = "Adjusted")
) 

man_beta_data_mean <- man_beta_data |>
  summarize(bias_est = mean(bias_est, na.rm = TRUE), .by = c("type", "cov"))

man_beta_data |>
  ggplot(aes(x = bias_est, color = type, linetype = type)) +
  geom_vline(xintercept = 0) +
  geom_vline(data = man_beta_data_mean, aes(xintercept = bias_est, color = type), linewidth = 1, linetype = 2) +
  geom_density(linewidth = 1) +
  facet_wrap(~factor(cov, c("Unadjusted", "Adjusted")), ncol = 1) +
  labs(
    title = "Bias in X1 beta coefficient for Y1 estimation under missing scenarios",
    subtitle = c("n = 1000; m = 5; miss_prop = 0.3; MAR dataset = 100; manual missingness"),
    x = "Bias (Beta_est - Beta_true)",
    y = "Density"
  ) +
  ms::scale_color_ms() +
  ms::theme_ms() +
  theme(legend.title = element_blank())
ggsave(
  "~/Downloads/man_beta_bias_simul.pdf",
  width = 7, height = 5, device = cairo_pdf
)
