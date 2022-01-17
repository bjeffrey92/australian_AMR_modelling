source("EDA.R")

load_inla_result <- function(species, ab) {
    inla_results <- readRDS(
        sprintf("results/inla/inla_result_%s_%s.rds", species, ab)
    )
    index <- as.integer(row.names(inla_results$validation_data))
    mean_prediction <- inla_results$result$summary.linear.predictor[
        index, "mean"
    ]

    validation_data <- inla_results$validation_data
    validation_data$inla_forecast <- mean_prediction
    validation_data$Forecasting_Horizon <- as.numeric(
        as.factor(
            validation_data$Year_Quarter
        )
    ) / 4
    validation_data$Location <- validation_data$Postcode
    headers <- c(
        "Location",
        "inla_forecast",
        "Forecasting_Horizon",
        "Perc.R"
    )
    return(validation_data[
        ,
        headers
    ])
}


load_ts_result <- function(species, ab) {
    ts_results <- read.csv(
        sprintf(
            "results/time_series/time_series_result_%s_%s.tsv",
            species,
            ab
        ),
        sep = "\t"
    )
    headers <- c("Location", "Forecasting_Horizon")
    headers <- c(
        headers,
        names(ts_results)[endsWith(names(ts_results), "cast")]
    )
    return(ts_results[, headers])
}


clip_predictions <- function(predictions, min = 0, max = 100) {
    predictions <- pmin(predictions, max)
    return(pmax(predictions, min))
}


load_results <- function(species, ab) {
    ts_result <- load_ts_result(species, ab)
    inla_result <- load_inla_result(species, ab)
    result <- merge(ts_result, inla_result)
    result <- pivot_longer(
        result,
        c(
            "naive_forecast_forecast",
            "ets_forecast_forecast",
            "arima_forecast_forecast",
            "inla_forecast"
        ),
        names_to = "forecast_type",
        values_to = "forecast"
    )
    result$forecast <- clip_predictions(result$forecast)
    return(result)
}


for (i in species_ab_combos) {
    species <- i[[1]]
    ab <- i[[2]]
    result <- load_results(species, ab)
    result$error <- abs(result$Perc.R - result$forecast)
}