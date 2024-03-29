library(dplyr)

source("data_loader.R")


exclude_infrequent <- function(df, threshold = 10) {
    df <- df %>%
        group_by(Postcode, Antibiotic) %>%
        filter(n() > threshold)
    return(df)
}


drop_locs_with_missing <- function(df) {
    n_time_points <- length(
        ts(
            start = min(df$Year_Quarter),
            end = max(df$Year_Quarter),
            frequency = 4
        )
    )
    inspect_time_series <- function(data) {
        if (
            nrow(data) == n_time_points
        ) {
            return(data)
        } else {
            data.frame(
                matrix(
                    ncol = ncol(data),
                    nrow = 0,
                    dimnames = list(NULL, names(data))
                )
            )
        }
    }
    tables <- lapply(
        group_by(df, Postcode, .add = TRUE) %>% group_split(),
        inspect_time_series
    )
    return(do.call("rbind", tables))
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


ecoli_abs <- c("AMPI", "AUGM", "CLEX", "TRIM")
staph_abs <- c("ERYT", "MUPI", "PENI", "TETR")

# for aggregated plots
all_ecoli_abs <- c("AMPI", "AUGM", "CLEX", "TRIM", "NORF")
all_staph_abs <- c("ERYT", "MUPI", "PENI", "TETR", "CLEX")

species_ab_combos <- c(
    lapply(ecoli_abs, function(x) c("Ecoli", x)),
    lapply(staph_abs, function(x) c("Staph", x))
)

paper_species_ab_combos <- c(
    lapply(all_ecoli_abs, function(x) c("Ecoli", x)),
    lapply(all_staph_abs, function(x) c("Staph", x))
)


validate_data <- function(df) {
    df <- drop_locs_with_missing(df)

    location_counts <- table(df$Postcode)
    if (length(location_counts) < 5 | min(location_counts) < 10) {
        return(NULL)
    }
    return(df)
}