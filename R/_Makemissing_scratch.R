# run top of 02_ampute_impute_explore.R first
expanded

patterns <- rbind(
    c(1,1,1,1,0,1,1,1,1,1,1,1),
    c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1),
    c(1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1),
    c(1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1)
)
patterns

f <- rep(1/nrow(patterns), nrow(patterns))

mar <- as_tibble(MAR(expanded, alpha = 0.5, pattern = patterns, f = f))
names(mar) <- names(expanded)
mnar <- as_tibble(MNAR(expanded, alpha = 0.5, pattern = patterns, f = f))
names(mnar) <- names(expanded)

glm(t2d ~ bmi, data = expanded, family = "binomial") |> summary()
glm(t2d ~ bmi, data = mar, family = "binomial") |> summary()
glm(t2d ~ bmi, data = mnar, family = "binomial") |> summary()
