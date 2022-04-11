library(forecast)

source("data_loader.R")
source("EDA.R")

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


check_seasonal <- function(species, ab) {
    df <- cache_load_data(species)
    df <- df[df$Antibiotic == ab, ]

    df <- validate_data(df)
    df <- do.call("rbind", by(df, df$Year_Quarter, aggregate_locations))

    time_series <- ts(df$Perc.R,
        start = min(df$Year_Quarter),
        frequency = 4
    )
    fit <- ets(time_series)
    return(c(species, ab, fit$method))
}

check_seasonal <- Vectorize(check_seasonal, c("ab"))