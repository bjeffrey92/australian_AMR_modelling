library(MyForecaster)

source("data_loader.R")
source("EDA.R")

forecasting_functions <- head(england_paper_forecasting_functions, 3)
do_forecasts_ <- Vectorize(do_forecasts, c("forecasting_horizon"))

run_forecasts <- function(data, organism, output_file, PI = 95) {
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
    write_output(t(output), output_file)
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


main <- function(species, ab, prediction_steps = 4, aggregate = FALSE) {
    df <- cache_load_data(species)
    df <- df[df$Antibiotic == ab, ]

    df <- validate_data(df)
    if (is.null(df)) {
        return(FALSE)
    }
    output_file <- sprintf("time_series_result_%s_%s.tsv", species, ab)
    if (aggregate) {
        df <- do.call("rbind", by(df, df$Year_Quarter, aggregate_locations))
    }
    output_file <- paste0("aggregated_", output_file)
    write_output_file_headers(output_file, forecasting_functions)
    for (i in group_split(df, Postcode)) {
        run_forecasts(i, species, output_file = output_file)
    }
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