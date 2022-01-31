library(MyINLA)

source("data_loader.R")
source("EDA.R")


beta_censoring_correction <- function(y_values) {
    min_value <- min(y_values[y_values != 0])
    y_values[y_values == 0] <- min_value / 1000
    y_values[y_values == 1] <- 1 - min_value / 1000
    return(y_values)
}


main <- function(species, ab, prediction_steps = 4, AR = 4, likelihood = "beta") {
    df <- cache_load_data(species)
    df <- df[df$Antibiotic == ab, ]

    df <- validate_data(df)
    if (is.null(df)) {
        return(FALSE)
    }

    df$prop.R <- df$Perc.R / 100

    # needs to start from 1
    df$time_point <- as.numeric(as.factor(df$Year_Quarter))

    x_col <- "long"
    y_col <- "lat"
    time_col <- "time_point"
    response_col <- "prop.R"
    location_col <- "Postcode"
    covariate_cols <- c()

    if (likelihood == "beta") {
        df[response_col] <- beta_censoring_correction(df[response_col])
    }

    if (prediction_steps > 0) {
        last_time_points <- tail(
            sort(unique(df[[time_col]])),
            prediction_steps
        )
        df <- as.data.frame(df) # convert to dataframe to retain the index
        data <- df[!df[[time_col]] %in% last_time_points, ]
        validation_data <- df[df[[time_col]] %in% last_time_points, ]

        mesh_and_spde <- build_mesh_and_spde(data, x_col, y_col, time_col,
            predict_forwards = prediction_steps, cutoff = 0.01
        )
    } else {
        mesh_and_spde <- build_mesh_and_spde(
            df,
            x_col,
            y_col,
            time_col,
            cutoff = 0.01
        )
        validation_data <- NULL
    }
    spde <- mesh_and_spde$spde

    stack_and_formula <- build_stack_and_formula(
        data, x_col, y_col, time_col,
        response_col, AR, covariate_cols,
        validation_data, prediction_steps,
        mesh_and_spde[["mesh2d"]],
        mesh_and_spde[["mesh1d"]]
    )
    stack.est <- stack_and_formula$stack.est
    stack <- stack_and_formula$stack


    formula <- as.formula(stack_and_formula[["formula_string"]])
    mod.mode <- inla(formula,
        data = inla.stack.data(stack.est, spde = spde),
        family = "gaussian",
        control.predictor = list(A = inla.stack.A(stack.est), compute = FALSE),
        control.compute = list(cpo = FALSE),
        keep = FALSE, verbose = TRUE,
        control.inla = list(reordering = "metis")
    )
    mod <- inla(formula,
        data = inla.stack.data(stack, spde = spde),
        family = "gaussian",
        control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
        control.compute = list(cpo = TRUE, dic = TRUE),
        control.mode = list(theta = mod.mode$mode$theta, restart = FALSE),
        keep = FALSE, verbose = TRUE,
        control.inla = list(reordering = "metis")
    )

    # result <- run_model(
    #     stack_and_formula[["stack"]],
    #     stack_and_formula[["formula_string"]],
    #     mesh_and_spde[["spde"]],
    #     likelihood = likelihood
    # )

    out <- list(
        "mesh_and_spde" = mesh_and_spde,
        "stack_and_formula" = stack_and_formula,
        "result" = result,
        "data" = data
    )
    if (prediction_steps > 0) {
        validation <- validate_model_fit(
            result,
            data,
            validation_data,
            response_col,
            time_col,
            location_col,
            stack_and_formula[["stack.est"]],
            prediction_steps
        )
        out[["validation"]] <- validation
        out[["validation_data"]] <- validation_data
    }
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