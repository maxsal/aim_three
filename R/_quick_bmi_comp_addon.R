d2 <- qread(glue("data/private/demo_comb_labs_ehr_20230322_EM_202.2.qs"))
pim2 <- qread(glue("data/private/mgi/20220822/MGI_PIM0X_20220822_EM_202.2.qs"))

comb <- left_join(
    d2 |>
        mutate(nhw = as.numeric(race_eth == "NHW")) |>
        select(id, age, female, nhw, smoke = smoke_ever, bmi = bmi_med, glucose = glucose_med),
    pim2 |>
        select(id, t2d = case),
    by = "id"
) |> select(-id)

mod1 <- glm(glucose ~ bmi, data = comb)

imp2 <- mice(comb |> select(all_of(c("age", "female", "smoke", "bmi", "glucose"))), m = 10, seed = 123, printFlag = FALSE)

pool(with(imp2, glm(glucose ~ bmi)))


comb_nhw <- comb |>
    filter(nhw == 1) |>
    mutate(
        bmi = scale(bmi)[, 1]
    )
mod2 <- glm(t2d ~ bmi, data = comb_nhw, family = "binomial")
imp3 <- mice(comb_nhw |> select(all_of(c("age", "female", "smoke", "bmi", "t2d"))), m = 10, seed = 123, printFlag = FALSE)
map_dbl(1:10,
    \(i) {
        glm(t2d ~ bmi, data = complete(imp3, action = i), family = "binomial")[["coefficients"]][["bmi"]]
    }) |> mean()

pool(with(imp3, glm(t2d ~ bmi, family = "binomial")))
