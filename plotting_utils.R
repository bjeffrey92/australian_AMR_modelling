source("EDA.R")

library(INLA)
library(stringr)

load_inla_result <- function(species,
                             ab,
                             likelihood,
                             percentage_correction = FALSE,
                             total_freq = FALSE,
                             get_PI = FALSE) {
    inla_results <- readRDS(
        sprintf(
            "results/inla/inla_result_%s_likelihood_%s_%s.rds",
            likelihood,
            species,
            ab
        )
    )


    validation_data <- inla_results$all_data[
        inla_results$all_data$train_test == "test",
    ]
    validation_data$inla_forecast <- validation_data$fitted_values
    if (percentage_correction) {
        validation_data$inla_forecast <- validation_data$inla_forecast * 100
    }

    validation_data$Forecasting_Horizon <- as.numeric(
        as.factor(
            validation_data$Year_Quarter
        )
    ) / 4
    validation_data$Location <- validation_data$Postcode
    headers <- c(
        "Location",
        "inla_forecast",
        "Forecasting_Horizon",
        "Perc.R"
    )
    if (total_freq) {
        headers <- c(headers, "totalFreq")
    }
    if (get_PI) {
        index_pred <- inla.stack.index(inla_results$stack, "pred")$data
        val_preds_upper <- inla_results$mod_mode$summary.fitted.val[
            index_pred, "0.975quant"
        ]
        val_preds_lower <- inla_results$mod_mode$summary.fitted.val[
            index_pred, "0.025quant"
        ]
        validation_data$fc_upper <- val_preds_upper
        validation_data$fc_lower <- val_preds_lower

        if (percentage_correction) {
            validation_data$fc_upper <- validation_data$fc_upper * 100
            validation_data$fc_lower <- validation_data$fc_lower * 100
        }

        headers <- c(headers, "fc_upper", "fc_lower")
    }
    return(validation_data[
        ,
        headers
    ])
}


load_ts_result <- function(species, ab, get_PI = FALSE) {
    ts_results <- read.csv(
        sprintf(
            "results/time_series/time_series_result_%s_%s.tsv",
            species,
            ab
        ),
        sep = "\t"
    )
    headers <- c("Location", "Forecasting_Horizon")
    forecast_headers <- names(ts_results)[endsWith(names(ts_results), "cast")]
    headers <- c(headers, forecast_headers)
    out <- ts_results[, headers]

    if (get_PI) {
        pi_widths <- ts_results[
            ,
            names(ts_results)[endsWith(names(ts_results), "PI_Width")]
        ]
        pis <- pi_widths / 2
        for (i in forecast_headers) {
            col_start <- str_split(i, "cast_forecast")[[1]][[1]]
            col <- paste0(col_start, "cast_PI_Width")
            out[[paste0(col_start, "cast_upper")]] <-
                ts_results[[i]] + pis[[col]]
            out[[paste0(col_start, "cast_lower")]] <-
                ts_results[[i]] - pis[[col]]
        }
    }

    return(out)
}


clip_predictions <- function(predictions, min = 0, max = 100) {
    predictions <- pmin(predictions, max)
    return(pmax(predictions, min))
}


apply_ensembling <- function(df, ensemble_model) {
    ensemble_rows <-
        df[
            Vectorize(startsWith, c("prefix"))
            (df$forecast_type, ensemble_model) %>% apply(1, any),
        ]
    ensemble_rows$forecast_type <- "ensemble_forecast"
    ensemble_rows$forecast <- mean(ensemble_rows$forecast)
    return(
        ensemble_rows[1, ]
    )
}


load_results <- function(species,
                         ab,
                         inla_likelihood,
                         percentage_correction = FALSE,
                         ensemble_model = c()) {
    ts_result <- load_ts_result(species, ab)
    inla_result <- load_inla_result(
        species,
        ab,
        inla_likelihood,
        percentage_correction
    )
    result <- merge(ts_result, inla_result)
    result <- pivot_longer(
        result,
        c(
            "naive_forecast_forecast",
            "ets_forecast_forecast",
            "arima_forecast_forecast",
            "inla_forecast"
        ),
        names_to = "forecast_type",
        values_to = "forecast"
    )
    result$forecast <- clip_predictions(result$forecast)

    if (length(ensemble_model) > 0) {
        ensemble_results <-
            result %>%
            group_split(Location, Forecasting_Horizon) %>%
            lapply(apply_ensembling, ensemble_model)
        ensemble_results <- do.call("rbind", ensemble_results)
        result <- rbind(as.data.frame(result), ensemble_results)
    }

    result$error <- abs(result$Perc.R - result$forecast)
    return(result)
}
load_results_vec <- Vectorize(load_results, c("ab"))


format_results <- function(result, ordered_forecasts = NULL) {
    rename_forecast_types <- function(f_types) {
        return(sapply(
            str_split(f_types, "_"),
            function(x) toupper(x[1])
        ))
    }
    result$forecast_type <- rename_forecast_types(result$forecast_type)
    if (!is.null(ordered_forecasts)) {
        result$forecast_type <-
            factor(
                result$forecast_type,
                levels = rename_forecast_types(ordered_forecasts),
                ordered = TRUE
            )
    }
    result$Forecasting_Horizon <- result$Forecasting_Horizon * 12
    result$Forecasting_Horizon <- as.factor(result$Forecasting_Horizon)
    names(result)[names(result) == "Forecasting_Horizon"] <-
        "Forecasting Horizon (months)"
    return(result)
}

order_forecast_accuracies <- function(result) {
    mean_forecast_error <- result %>%
        group_split(forecast_type) %>%
        sapply(function(x) c(mean(x$error), x$forecast_type[[1]])) %>%
        t()
    return(mean_forecast_error[order(mean_forecast_error[, 1]), 2])
}


load_aggregated_results <- function(species,
                                    ab,
                                    models = c("inla", "arima", "ets")) {
    df <- read.csv(
        sprintf(
            "results/aggregated_locations_results/aggregated_time_series_result_%s_%s.tsv",
            species,
            ab
        ),
        sep = "\t"
    )
    parse_model_fc <- function(model) {
        fc <- df[, paste0(model, "_forecast_forecast")]
        fc_upper <- sapply(
            fc + df[, paste0(model, "_forecast_PI_Width")] / 2,
            function(x) min(x, 100)
        )
        fc_lower <- sapply(
            fc - df[, paste0(model, "_forecast_PI_Width")] / 2,
            function(x) max(x, 0)
        )
        return(data.frame(fc, fc_upper, fc_lower, model))
    }
    data <- lapply(models[models != "inla"], parse_model_fc)

    if ("inla" %in% models) {
        inla_result <- load_inla_result(
            species,
            ab,
            "beta",
            percentage_correction = FALSE,
            total_freq = TRUE,
            get_PI = TRUE
        )

        aggregate_col <- function(totalFreq, col) {
            return(sum(col * totalFreq) / sum(totalFreq) * 100)
        }
        df <- group_by(inla_result, Forecasting_Horizon) %>%
            summarise(
                fc_upper = aggregate_col(totalFreq, fc_upper),
                fc_lower = aggregate_col(totalFreq, fc_lower),
                fc = aggregate_col(totalFreq, inla_forecast)
            )
        df$model <- "inla"
        df <- as.data.frame(df)
        data[[length(data) + 1]] <- df[, c("fc", "fc_upper", "fc_lower", "model")]
    }
    return(data)
}


tabulate_spatial_effect <- function(species, abs_list) {
    spatial_effects <- lapply(abs_list, load_spatial_field, species)
    spatial_effects <- do.call("rbind", spatial_effects)
    spatial_effects <- spatial_effects[order(row.names(spatial_effects)), ]
    range_and_var <- extract_range_and_var(spatial_effects)
    range_and_var <-
        cbind(
            head(range_and_var, length(abs_list)),
            tail(range_and_var, length(abs_list))
        )
    df <- as.data.frame(
        range_and_var,
        row.names = abs_list,
    )
    names(df) <- c("range", "variance")
    df$variance <- str_replace(df$variance, "\\*", "")
    write.csv(df, sprintf("%s_spatial_effects.csv", species))
}


add_ab_and_species <- function(results, abs_list, species) {
    row_counts <- sapply(results, nrow)
    results <- do.call("rbind", results)
    results$species <- species
    results$ab <- rep(abs_list, row_counts)
    return(results[results$forecast_type == "inla_forecast", ])
}