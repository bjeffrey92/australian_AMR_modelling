library(dplyr)
library(memoise)
library(Microsoft365R)
library(zoo)

format_data <- function(df) {
    df <- dplyr::tibble(df)
    df$Year_Quarter <- zoo::as.yearqtr(df$Year_Quarter)
    df <- dplyr::arrange(df, df$Year_Quarter)
    df <- df[order(df$Year_Quarter), ]
    return(df)
}


postcode_lat_longs <- function() {
    lat_long <- read.csv(
        "https://www.matthewproctor.com/Content/postcodes/australian_postcodes.csv"
    )
    return(tibble(lat_long[, c("postcode", "lat", "long")]))
}


load_data <- function(species) {
    ob <- Microsoft365R::get_business_onedrive(auth_type = "device_code")
    fname <- sprintf(
        "ab_time_series_stuff/Australian_Data/%s_greater10.csv",
        species
    )
    save_location <- tempfile(pattern = "", fileext = ".csv")
    ob$download_file(fname, save_location)
    df <- read.csv(save_location)
    unlink(save_location)

    lat_long <- postcode_lat_longs()
    lat_long <- lat_long[lat_long$postcode %in% df$Postcode, ]
    lat_long <- group_by(lat_long, postcode) %>%
        mutate(lat = mean(lat), long = mean(long)) %>%
        unique()
    df <- merge(df, lat_long, by.x = "Postcode", by.y = "postcode") %>%
        tibble()

    return(format_data(df))
}

cache_load_data <- memoise(load_data)