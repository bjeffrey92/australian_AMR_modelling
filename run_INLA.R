source("INLA_driver.R")
source("data_loader.R")
source("EDA.R")


beta_censoring_correction <- function(y_values) {
    min_value <- min(y_values[y_values != 0])
    y_values[y_values == 0] <- min_value / 1000
    y_values[y_values == 1] <- 1 - min_value / 1000
    return(y_values)
}


import_and_format_data <- function(species, ab, prediction_steps, likelihood) {
    df <- cache_load_data(species)
    df <- df[df$Antibiotic == ab, ]

    df <- validate_data(df)
    if (is.null(df)) {
        return(NULL)
    }

    # needs to start from 1
    df$time_point <- as.numeric(as.factor(df$Year_Quarter))
    df$prop.R <- df$Perc.R / 100
    if (likelihood == "beta") {
        df$prop.R <- beta_censoring_correction(df$prop.R)
    }

    last_time_points <- tail(
        sort(unique(df[[time_col]])),
        prediction_steps
    )
    df <- as.data.frame(df) # convert to dataframe to retain the index
    data <- df[!df$time_point %in% last_time_points, ]
    validation_data <- df[df$time_point %in% last_time_points, ]

    return(list(data, validation_data))
}


main <- function(species, ab, prediction_steps = 4,
                 AR = 4, likelihood = "beta") {
    x_col <- "long"
    y_col <- "lat"
    time_col <- "time_point"
    response_col <- "prop.R"
    location_col <- "Postcode"

    all_data <- import_and_format_data(
        species,
        ab,
        prediction_steps,
        likelihood
    )
    if (is.null(all_data)) {
        return(FALSE)
    }
    data <- all_data[[1]]
    validation_data <- all_data[[2]]

    spde <- build_mesh_and_spde(data, x_col, y_col, time_col)

    stacks <- build_stack(
        data, validation_data, x_col,
        y_col, time_col, response_col, spde
    )
    stack_est <- stacks[[1]]
    stack_pred <- stacks[[2]]
    stack <- inla.stack(stack_est, stack_pred)

    formula_string <-
        sprintf(
            "outcome ~ -1 + Intercept + f(field, model = spde, group = field.group, control.group = list(model = 'ar', order=%s))",
            AR
        )
    formula <- as.formula(formula_string)

    mod_mode <- inla(formula,
        data = inla.stack.data(stack_est, spde = spde),
        family = likelihood,
        control.predictor = list(A = inla.stack.A(stack_est), compute = FALSE),
        control.compute = list(cpo = FALSE),
        keep = FALSE, verbose = TRUE,
        control.inla = list(reordering = "metis")
    )
    mod <- inla(formula,
        data = inla.stack.data(stack, spde = spde),
        family = likelihood,
        control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
        control.compute = list(cpo = TRUE, dic = TRUE),
        control.mode = list(theta = mod_mode$mode$theta, restart = FALSE),
        keep = FALSE, verbose = TRUE,
        control.inla = list(reordering = "metis")
    )

    index_est <- inla.stack.index(stack, "est")$data
    index_pred <- inla.stack.index(stack, "pred")$data
    est_predictions <- mod$summary.linear.predictor[index_est, "mean"]
    val_predictions <- mod$summary.linear.predictor[index_pred, "mean"]


    # validation <- validate_model_fit(
    #     result,
    #     data,
    #     validation_data,
    #     response_col,
    #     time_col,
    #     location_col,
    #     stack_and_formula[["stack.est"]],
    #     prediction_steps
    # )
    out <- list(
        "mesh_and_spde" = mesh_and_spde,
        "stack_and_formula" = stack_and_formula,
        "result" = result,
        "data" = data,
        "validation" = validation,
        "validation_data" <- validation_data
    )

    saveRDS(
        out,
        sprintf("inla_result_%s_likelihood_%s_%s.rds", likelihood, species, ab)
    )
    return(TRUE)
}

for (i in species_ab_combos) {
    success <- main(i[[1]], i[[2]])
    if (success) {
        sprintf("%s completed", i)
    } else {
        sprintf("insufficient data for %s", i)
    }
}