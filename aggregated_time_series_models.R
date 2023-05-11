library(MyForecaster)

source("data_loader.R")
source("EDA.R")

forecasting_functions <- head(england_paper_forecasting_functions, 3)
do_forecasts_ <- Vectorize(do_forecasts, c("forecasting_horizon"))

run_forecasts <- function(data, organism, PI = 95) {
    ab <- data$Antibiotic[[1]]
    location <- data$Postcode[[1]]

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


main <- function() {
    all_forecasts <- lapply(paper_species_ab_combos, gather_all_forecasts)
}