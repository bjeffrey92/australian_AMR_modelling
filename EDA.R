library(dplyr)

source("data_loader.R")


exlcude_infrequent <- function(df, threshold = 10) {
    df <- df %>%
        group_by(Postcode, Antibiotic) %>%
        filter(n() > threshold)
    return(df)
}


aggregate_freq <- function(quart_data) {
    return(
        data.frame(
            totalFreq = sum(quart_data$R) / sum(quart_data$totalFreq) * 100,
            Year_Quarter = quart_data$Year_Quarter[[1]],
            Antibiotic = quart_data$Antibiotic[[1]]
        )
    )
}


plot_time_series <- function(species) {
    df <- load_data(species)
    aggregated_data <- lapply(
        group_split(df, Antibiotic, Year_Quarter),
        aggregate_freq
    )
    aggregated_data <- do.call(rbind, aggregated_data)
    for (data in group_split(aggregated_data, Antibiotic)) {
        ab <- data$Antibiotic[[1]]
        png(sprintf("%s_%s_aggregated_locations.png", species, ab))
        plot(
            data$Year_Quarter,
            data$totalFreq,
            type = "l",
            main = sprintf("percent resistance %s %s", ab, species)
        )
        dev.off()
    }
}
