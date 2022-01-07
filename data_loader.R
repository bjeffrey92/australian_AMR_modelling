library(Microsoft365R)

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
    return(df)
}
