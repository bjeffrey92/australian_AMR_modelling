library(MyForecaster)
library(purrr)

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
    ets_fc <- unlist(output["ets_forecast_forecast", ])
    arima_fc <- unlist(output["arima_forecast_forecast", ])

    fc_list <- list(naive_fc = naive_fc, ets_fc = ets_fc, arima_fc = arima_fc)
    fc_list <- lapply(fc_list, ts, start = max(data$Year_Quarter) - 0.75, frequency = 4)

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


ts_to_dataframe <- function(data) {
    df <- data.frame(Y = as.matrix(data), date = time(data))
    return(df)
}


format_forecast_data <- function(forecast) {
    df <- data.frame(
        time_series = as.matrix(forecast$time_series),
        date = time(forecast$time_series)
    )
    forecasts <- lapply(forecast$fc_list, ts_to_dataframe) %>%
        reduce(left_join, by = "date")
    names(forecasts)[names(forecasts) == "Y.x"] <- "Naive"
    names(forecasts)[names(forecasts) == "Y.y"] <- "ETS"
    names(forecasts)[names(forecasts) == "Y"] <- "ARIMA"

    df <- left_join(df, forecasts, by = "date")
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
}