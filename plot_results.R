source("EDA.R")

library(cowplot)
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
    validation_data <- inla_results$all_data[
        inla_results$all_data$train_test == "test",
    ]
    validation_data$inla_forecast <- validation_data$fitted_values
    if (percentage_correction) {
        validation_data$inla_forecast <- validation_data$inla_forecast * 100
    }

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


load_results <- function(species, ab, inla_likelihood, percentage_correction = FALSE) {
    ts_result <- load_ts_result(species, ab)
    inla_result <- load_inla_result(
        species,
        ab,
        inla_likelihood,
        percentage_correction
    )
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


format_results <- function(result, ordered_forecasts = NULL) {
    rename_forecast_types <- function(f_types) {
        return(sapply(
            str_split(f_types, "_"),
            function(x) paste(head(x, 2), collapse = " ")
        ))
    }
    result$forecast_type <- rename_forecast_types(result$forecast_type)
    if (!is.null(ordered_forecasts)) {
        result$forecast_type <-
            factor(
                result$forecast_type,
                levels = rename_forecast_types(ordered_forecasts),
                ordered = TRUE
            )
    }
    result$Forecasting_Horizon <- as.factor(result$Forecasting_Horizon)
    return(result)
}


plot_result <- function(result, species, ab,
                        ordered_forecasts = NULL, save = TRUE) {
    result <- format_results(result, ordered_forecasts)
    p <- ggplot(
        result,
        aes(x = forecast_type, y = error, fill = forecast_type)
    ) +
        geom_boxplot() +
        theme(
            # text = element_text(size = 22),
            axis.text.x = element_text(angle = 60, hjust = 1)
        ) +
        theme_minimal()
    if (save) {
        ggsave(sprintf("%s_%s_forecast_accuracy.png", species, ab),
            plot = p,
        )
    } else {
        return(p)
    }
}


order_forecast_accuracies <- function(result) {
    mean_forecast_error <- result %>%
        group_split(forecast_type) %>%
        sapply(function(x) c(mean(x$error), x$forecast_type[[1]])) %>%
        t()
    return(mean_forecast_error[order(mean_forecast_error[, 1]), 2])
}


plot_time_series <- function(df, ab, result) {
    df <- df[df$Antibiotic == ab, ]
    df <- df[order(df$Year_Quarter), ]
    p <- ggplot(
        df,
        aes(x = Year_Quarter, y = totalFreq)
    ) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        xlab("Date") +
        ylab("Percent Resistant")
    return(p)
}


time_series_accuracy_plotter <- function(species, ab) {
    result <- load_results(
        species,
        ab,
        inla_likelihood = "beta",
        percentage_correction = TRUE
    )
    result$error <- abs(result$Perc.R - result$forecast)

    ordered_forecasts <- order_forecast_accuracies(result)
    accuracy_plot <- plot_result(
        result,
        species,
        ab,
        ordered_forecasts = ordered_forecasts,
        save = FALSE
    )

    time_series_data <- cache_load_data(species)
    aggregated_data <- lapply(
        group_split(time_series_data, Antibiotic, Year_Quarter),
        aggregate_freq
    )
    aggregated_data <- do.call(rbind, aggregated_data)
    time_series_plot <- plot_time_series(aggregated_data)
}


panel_plots <- function(species, abs) {

}



for (i in species_ab_combos) {
    species <- i[[1]]
    ab <- i[[2]]
    tryCatch(
        expr = {
            result <- load_results(
                species,
                ab,
                inla_likelihood = "beta",
                percentage_correction = TRUE
            )
            result$error <- abs(result$Perc.R - result$forecast)
            plot_result(result, species, ab)
        },
        error = function(e) {
            print(e)
            sprintf(
                "Failed for %s %s", species, ab
            )
        }
    )
}