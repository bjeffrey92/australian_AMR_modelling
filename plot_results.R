source("EDA.R")

library(cowplot)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggpubr)


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
            function(x) toupper(x[1])
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
        aes(
            x = forecast_type,
            y = error,
            # group = interaction(forecast_type, Forecasting_Horizon),
            fill = Forecasting_Horizon
        )
    ) +
        geom_boxplot() +
        theme_minimal() +
        theme(legend.position = "bottom") +
        xlab("Forecast Type") +
        ylab("Forecast Error")
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


plot_time_series <- function(df, ab, forecast_horizon) {
    df <- df[df$Antibiotic == ab, ]
    df <- df[order(df$Year_Quarter), ]
    p <- ggplot(
        df,
        aes(x = Year_Quarter, y = totalFreq)
    ) +
        geom_line() +
        geom_point() +
        geom_vline(
            xintercept = max(df$Year_Quarter) - forecast_horizon,
            linetype = "dashed"
        ) +
        theme_minimal() +
        xlab("Date") +
        ylab("Total % Resistant")
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
    forecast_horizon <- max(result$Forecasting_Horizon)
    time_series_plot <- plot_time_series(aggregated_data, ab, forecast_horizon)
    return(list(time_series_plot, accuracy_plot))
}
time_series_accuracy_plotter <- Vectorize(time_series_accuracy_plotter, c("ab"))


panel_plots <- function(species, abs_list, save_legend = FALSE) {
    plots <- time_series_accuracy_plotter(species, abs_list)

    if (save_legend) {
        legend <- cowplot::get_legend(plots[[2]]) %>% as_ggplot()
        ggsave("legend.png", plot = legend)
    }
    plots <- lapply(plots, function(x) x + theme(legend.position = "none"))

    plots <- c(
        lapply(plots[1:6], function(x) x + xlab("")), plots[7:8]
    )
    if (species == "Ecoli") {
        species_ <- "E. coli"
    } else if (species == "Staph") {
        species_ <- "S. aureus"
    }

    panel_plot <- cowplot::plot_grid(
        plotlist = plots,
        ncol = 2,
        labels = c(
            rbind(paste(sprintf("%s -", species_), abs_list), rep("", 4))
        ),
        label_x = 0.1,
        label_y = 1.1,
        label_size = 9
    ) +
        theme(plot.margin = margin(20))

    ggsave(sprintf("%s_forecast_panel_plot.png", species),
        plot = panel_plot
    )
}


load_spatial_field <- function(ab, species, likelihood = "beta") {
    inla_results <- readRDS(
        sprintf(
            "results/inla/inla_result_%s_likelihood_%s_%s.rds",
            likelihood,
            species,
            ab
        )
    )
    hyperparams <- inla_results$mod_mode$summary.hyperpar
    theta <- hyperparams[2:3, 3:5]
    theta$significant <- apply(theta > 0, 1, function(x) var(x) == 0)
    theta <- theta %>%
        mutate(asterisk = case_when(
            significant ~ "*",
            !significant ~ ""
        ))
    theta[1:3] <- theta[1:3] %>% round(2)
    return(
        sprintf(
            "%s (%s)%s",
            theta[[2]],
            paste(theta[[1]], theta[[3]], sep = "-"),
            theta[[5]]
        )
    )
}


tabulate_spatial_effect <- function(species, abs_list) {
    spatial_effects <- lapply(abs_list, load_spatial_field, species)
    names(spatial_effects) <- abs_list
    df <- as.data.frame(do.call(rbind, spatial_effects))
    names(df) <- c("Theta_1", "Theta_2")
    write.csv(df, sprintf("%s_spatial_effects.csv", species))
}


panel_plots("Ecoli", ecoli_abs)
panel_plots("Staph", staph_abs)
tabulate_spatial_effect("Ecoli", ecoli_abs)
tabulate_spatial_effect("Staph", staph_abs)
