library(MyForecaster)
library(purrr)
library(tidyverse)
library(ggplot2)

source("data_loader.R")
source("EDA.R")

forecasting_functions <- head(england_paper_forecasting_functions, 3)
do_forecasts_ <- Vectorize(do_forecasts, c("forecasting_horizon"))

run_forecasts <- function(data, organism, PI = 95) {
    ab <- data$Antibiotic[[1]]
    location <- data$Postcode[[1]]

    # find position where the gap is more than 1 quarter if exists
    gaps <- diff(data$Year_Quarter) > 0.25
    if (sum(gaps) > 0) {
        jump_index <- which(gaps, TRUE)
        data <- head(data, jump_index)
    }

    time_series <- ts(data$Perc.R,
        start = min(data$Year_Quarter),
        frequency = 4
    )
    horizons <- 1:4 / 4
    output <- do_forecasts_(
        forecasting_functions,
        time_series,
        forecasting_horizon = horizons,
        max_forecasting_horizon = 4 / 4,
        location = location,
        ab = ab,
        PI = PI,
        organism = organism,
        frequency = 4
    )

    naive_fc <- unlist(output["naive_forecast_forecast", ])
    naive_fc_PI <- unlist(output["naive_forecast_PI_width", ]) / 2
    naive_fc_lower <- naive_fc - naive_fc_PI
    naive_fc_upper <- naive_fc + naive_fc_PI

    ets_fc <- unlist(output["ets_forecast_forecast", ])
    ets_fc_PI <- unlist(output["ets_forecast_PI_width", ]) / 2
    ets_fc_lower <- ets_fc - ets_fc_PI
    ets_fc_upper <- ets_fc + ets_fc_PI

    arima_fc <- unlist(output["arima_forecast_forecast", ])
    arima_fc_PI <- unlist(output["arima_forecast_PI_width", ]) / 2
    arima_fc_lower <- arima_fc - arima_fc_PI
    arima_fc_upper <- arima_fc + arima_fc_PI


    fc_list <- list(
        naive_fc = naive_fc,
        naive_fc_lower = naive_fc_lower,
        naive_fc_upper = naive_fc_upper,
        ets_fc = ets_fc,
        ets_fc_lower = ets_fc_lower,
        ets_fc_upper = ets_fc_upper,
        arima_fc = arima_fc,
        arima_fc_lower = arima_fc_lower,
        arima_fc_upper = arima_fc_upper
    )

    final_observation <- head(tail(time_series, 5), 1)
    fc_list <- lapply(fc_list, function(x) c(final_observation, x))

    fc_list <- lapply(
        fc_list, ts,
        start = max(data$Year_Quarter) - 1, frequency = 4
    )

    output_list <- list(time_series = time_series, fc_list = fc_list)
    return(output_list)
}


aggregate_locations <- function(df) {
    row <- df[1, ]
    row$Postcode <- -999
    row$lat <- 0
    row$long <- 0
    row$R <- sum(df$R)
    row$S <- sum(df$S)
    row$totalFreq <- sum(df$totalFreq)
    row$Perc.R <- row$R / row$totalFreq * 100
    return(row)
}


gather_all_forecasts <- function(species_ab) {
    species <- species_ab[[1]]
    ab <- species_ab[[2]]

    df <- cache_load_data(species)
    df <- df[df$Antibiotic == ab, ]

    df <- do.call("rbind", by(df, df$Year_Quarter, aggregate_locations))

    forecasts <- run_forecasts(df, species)
    forecasts$species <- species
    forecasts$ab <- ab

    return(forecasts)
}


format_forecast_data <- function(forecast) {
    df <- data.frame(
        time_series = as.matrix(forecast$time_series),
        date = time(forecast$time_series)
    )

    forecasts <- as.data.frame(do.call(cbind, forecast$fc_list))
    forecasts$date <- time(forecast$fc_list[[1]])

    df <- left_join(df, forecasts, by = "date")
}


panel_plot <- function(df, species) {
    names(df)[names(df) == "ets_fc"] <- "ETS"
    names(df)[names(df) == "arima_fc"] <- "ARIMA"

    colors <- c("ETS" = "red", "ARIMA" = "blue")

    p <- ggplot(data = df) +
        geom_line(aes(x = date, y = time_series)) +
        geom_point(aes(x = date, y = time_series)) +
        geom_line(aes(x = date, y = ETS, color = "ETS")) +
        geom_ribbon(
            aes(x = date, ymin = ets_fc_lower, ymax = ets_fc_upper),
            alpha = 0.1,
            fill = "red"
        ) +
        geom_line(aes(x = date, y = ARIMA, color = "ARIMA")) +
        geom_ribbon(
            aes(x = date, ymin = arima_fc_lower, ymax = arima_fc_upper),
            alpha = 0.1,
            fill = "blue"
        ) +
        facet_wrap(vars(ab), scales = "free") +
        theme_minimal() +
        ylab("Total % Resistant") +
        xlab("Date") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = c(0.8, 0.2)
        ) +
        scale_color_manual(values = colors) +
        labs(color = "Forecasting Model")

    # hard code this for now
    if (species == "Ecoli") {
        p <- p +
            geom_vline(
                data = filter(df, ab == "NORF"), aes(xintercept = 2016.5),
                linetype = "dashed"
            ) +
            geom_vline(
                data = filter(df, ab == "TRIM"), aes(xintercept = 2017.75),
                linetype = "dashed"
            ) +
            geom_vline(
                data = filter(df, ab == "AMPI"), aes(xintercept = 2017.75),
                linetype = "dashed"
            ) +
            geom_vline(
                data = filter(df, ab == "AUGM"), aes(xintercept = 2017.75),
                linetype = "dashed"
            ) +
            geom_vline(
                data = filter(df, ab == "CLEX"), aes(xintercept = 2017.75),
                linetype = "dashed"
            )
    } else {
        p <- p + geom_vline(
            data = df, aes(xintercept = 2017.75), linetype = "dashed"
        )
    }

    return(p)
}


main <- function() {
    names(paper_species_ab_combos) <- lapply(
        paper_species_ab_combos,
        paste0,
        collapse = "_"
    )
    all_forecasts <- lapply(paper_species_ab_combos, gather_all_forecasts)
    formatted_data <- lapply(all_forecasts, format_forecast_data)
    formatted_data <- Map(
        cbind,
        formatted_data,
        species_ab = names(formatted_data)
    )
    data <- do.call(rbind, formatted_data)
    data <- separate(
        data = data, col = species_ab, into = c("species", "ab"), sep = "_"
    )

    ecoli_plot <- panel_plot(data[data$species == "Ecoli", ], "Ecoli")
    staph_plot <- panel_plot(data[data$species == "Staph", ], "Staph")
}