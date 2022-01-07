library(dplyr)
library(Microsoft365R)
library(zoo)

format_data <- function(df) {
    df <- dplyr::tibble(df)
    df$Year_Quarter <- zoo::as.yearqtr(df$Year_Quarter)
    df <- dplyr::arrange(df, df$Year_Quarter)
    return(df)
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
    df <- format_data(df)
    return(df)
}
