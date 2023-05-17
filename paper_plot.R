library(tidyverse)

source("data_loader.R")


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

load_data <- function(species) {
    df <- cache_load_data(species)
    df <- do.call(
        "rbind",
        by(
            df, list(df$Year_Quarter, df$Antibiotic),
            aggregate_locations
        )
    )
    return(df)
}

make_summary_plot <- function(df) {
    names(df)[names(df) == "R"] <- "Resistant"
    names(df)[names(df) == "S"] <- "Sensitive"

    df2 <- gather(df, R_or_S, value, Resistant:Sensitive)
    p <- ggplot(data = df2) +
        geom_bar(aes(x = Year_Quarter, y = value, fill = R_or_S),
            position = "stack", stat = "identity"
        ) +
        facet_wrap(vars(Antibiotic)) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = c(0.8, 0.2)
        ) +
        xlab("Date") +
        ylab("# Isolates Sampled") +
        labs(fill = "Resistance Status")
    return(p)
}


main <- function() {
    ecoli_data <- load_data("Ecoli")
    staph_data <- load_data("Staph")
    staph_data <- staph_data[
        staph_data$Antibiotic %in% c("ERYT", "MUPI", "PENI", "TETR", "CLEX"), 
    ]

    ecoli_plot <- make_summary_plot(ecoli_data)
    staph_plot <- make_summary_plot(staph_data)
}