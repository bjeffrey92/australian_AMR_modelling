library(MyINLA)

source("data_loader.R")
source("EDA.R")

main <- function(species, ab, prediction_steps = 4, AR = 4) {
    df <- cache_load_data(species)
    df <- df[df$Antibiotic == ab, ]
    df <- drop_locs_with_missing(df)

    # needs to start from 1
    df$time_point <- as.numeric(as.factor(df$Year_Quarter))

    x_col <- "long"
    y_col <- "lat"
    time_col <- "time_point"
    response_col <- "totalFreq"
    location_col <- "Postcode"
    covariate_cols <- c()

    if (prediction_steps > 0) {
        last_time_points <- tail(
            sort(unique(data[[time_col]])),
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

    stack_and_formula <- build_stack_and_formula(
        data, x_col, y_col, time_col,
        response_col, AR, covariate_cols,
        validation_data, prediction_steps,
        mesh_and_spde[["mesh2d"]],
        mesh_and_spde[["mesh1d"]]
    )

    result <- run_model(
        stack_and_formula[["stack"]],
        stack_and_formula[["formula_string"]],
        mesh_and_spde[["spde"]]
    )

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
    saveRDS(out, sprintf("inla_result_%s_%s", species, ab))
}

species_ab_combos <- list(
    c("Ecoli", "AMPI"),
    c("Ecoli", "AUGM"),
    c("Ecoli", "CLEX"),
    c("Ecoli", "NORF"),
    c("Ecoli", "TRIM"),
    c("Staph", "CFOX"),
    c("Staph", "CLEX"),
    c("Staph", "COTR"),
    c("Staph", "ERYT"),
    c("Staph", "METH"),
    c("Staph", "MUPI"),
    c("Staph", "PENI"),
    c("Staph", "TETR")
)

for (i in species_ab_combos) {
    main(i[[1]], i[[2]])
}