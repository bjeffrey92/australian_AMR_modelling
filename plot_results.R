source("EDA.R")

library(ggplot2)
library(tidyr)
library(stringr)

load_inla_result <- function(species, ab, likelihood, percentage_correction = FALSE) {
    inla_results <- readRDS(
        sprintf(
            "results/inla/inla_result_%s_likelihood_%s_%s.rds",
            likelihood,
            species,
            ab
        )
    )
    index <- as.integer(row.names(inla_results$validation_data))
    mean_prediction <- inla_results$result$summary.linear.predictor[
        index, "mean"
    ]
    if (percentage_correction) {
        mean_prediction <- mean_prediction * 100
    }

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


load_results <- function(species, ab, inla_likelihood) {
    ts_result <- load_ts_result(species, ab)
    inla_result <- load_inla_result(species, ab, inla_likelihood)
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


format_results <- function(result) {
    result$forecast_type <- sapply(
        str_split(result$forecast_type, "_"),
        function(x) paste(head(x, 2), collapse = " ")
    )
    result$Forecasting_Horizon <- as.factor(result$Forecasting_Horizon)
    return(result)
}


plot_result <- function(result, species, ab) {
    result <- format_results(result)
    png(
        sprintf("%s_%s_forecast_accuracy.png", species, ab),
        width = 1080,
        height = 1080
    )
    ggplot(
        result,
        aes(x = forecast_type, y = error, fill = Forecasting_Horizon)
    ) +
        geom_boxplot() +
        theme(
            text = element_text(size = 22),
            axis.text.x = element_text(angle = 60, hjust = 1)
        )
    dev.off()
}

for (i in species_ab_combos) {
    species <- i[[1]]
    ab <- i[[2]]
    result <- load_results(species, ab, inla_likelihood = "gaussian")
    result$error <- abs(result$Perc.R - result$forecast)
    plot_result(result, species, ab)
}