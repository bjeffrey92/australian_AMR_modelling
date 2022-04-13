source("EDA.R")

library(cowplot)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggpubr)
library(GGally)
library(viridis)


plot_result <- function(result, species, ab,
                        ordered_forecasts = NULL, save = TRUE) {
    result <- format_results(result, ordered_forecasts)
    p <- ggplot(
        result,
        aes_string(
            x = "forecast_type",
            y = "error",
            fill = "`Forecasting Horizon (months)`"
        )
    ) +
        geom_boxplot() +
        theme_minimal() +
        theme(legend.position = "bottom") +
        xlab("Forecast Type") +
        ylab("Forecast Error")
    if (save) {
        ggsave(sprintf("%s_%s_forecast_accuracy.png", species, ab),
            plot = p,
        )
    } else {
        return(p)
    }
}


plot_time_series <- function(df, ab, forecast_horizon, species, plot_fits = TRUE) {
    df <- df[df$Antibiotic == ab, ]
    df <- df[order(df$Year_Quarter), ]
    p <- ggplot() +
        geom_line(data = df,
            aes(x = Year_Quarter, y = totalFreq)) +
        geom_point(data = df,
            aes(x = Year_Quarter, y = totalFreq)) +
        geom_vline(
            xintercept = max(df$Year_Quarter) - forecast_horizon,
            linetype = "dashed"
        ) +
        theme_minimal() +
        xlab("Date") +
        ylab("Total % Resistant")
    if (plot_fits) {
        forecasts <- load_aggregated_results(species, ab)
        dates <- as.yearqtr(
            seq(df$Year_Quarter[[nrow(df) - 4]],
                df$Year_Quarter[[nrow(df) - 4]] + 1, 0.25
            )
        )
        last_value <- rep(df$totalFreq[[nrow(df) - forecast_horizon * 4]], 4)
        for (model_fc in forecasts) {
            if (model_fc$model[[1]] == "ets") {
                col <- "red"
            } else if (model_fc$model[[1]] == "arima") {
               col <- "blue"
            }
            model_fc <- rbind(last_value, model_fc)
            model_fc$Year_Quarter <- dates
            p <- p +
                geom_line(data = model_fc,
                    aes(x = Year_Quarter, y = fc), col = col) +
                geom_ribbon(data = model_fc,
                    aes(x = Year_Quarter, ymin = fc_lower, ymax = fc_upper),
                    alpha = 0.1,
                    fill = col)
        }
    }
    return(p)
}


time_series_accuracy_plotter <- function(species, ab) {
    result <- load_results(
        species,
        ab,
        inla_likelihood = "beta",
        percentage_correction = TRUE
    )

    ordered_forecasts <- order_forecast_accuracies(result)
    accuracy_plot <- plot_result(
        result,
        species,
        ab,
        ordered_forecasts = ordered_forecasts,
        save = FALSE
    )

    time_series_data <- cache_load_data(species)
    aggregated_data <- lapply(
        group_split(time_series_data, Antibiotic, Year_Quarter),
        aggregate_freq
    )
    aggregated_data <- do.call(rbind, aggregated_data)
    forecast_horizon <- max(result$Forecasting_Horizon)
    time_series_plot <- plot_time_series(
        aggregated_data,
        ab,
        forecast_horizon,
        species
    )
    return(list(time_series_plot, accuracy_plot))
}
time_series_accuracy_plotter <- Vectorize(time_series_accuracy_plotter, c("ab"))


panel_plots <- function(species, abs_list, save_legend = FALSE) {
    plots <- time_series_accuracy_plotter(species, abs_list)

    if (save_legend) {
        legend <- cowplot::get_legend(plots[[2]]) %>% as_ggplot()
        ggsave("legend.png", plot = legend)
        return(NULL)
    }
    plots <- lapply(plots, function(x) x + theme(legend.position = "none"))

    plots <- c(
        lapply(plots[1:6], function(x) x + xlab("")), plots[7:8]
    )
    if (species == "Ecoli") {
        species_ <- "E. coli"
    } else if (species == "Staph") {
        species_ <- "S. aureus"
    }

    panel_plot <- cowplot::plot_grid(
        plotlist = plots,
        ncol = 2,
        labels = c(
            rbind(paste(sprintf("%s -", species_), abs_list), rep("", 4))
        ),
        label_x = 0.1,
        label_y = 1.1,
        label_size = 9
    ) +
        theme(plot.margin = margin(20))

    ggsave(sprintf("%s_forecast_panel_plot.png", species),
        plot = panel_plot, bg = "white"
    )
}


load_spatial_field <- function(ab, species, likelihood = "beta") {
    inla_results <- readRDS(
        sprintf(
            "results/inla/inla_result_%s_likelihood_%s_%s.rds",
            likelihood,
            species,
            ab
        )
    )
    result_field <- inla.spde.result(
        inla_results$mod_mode,
        "field",
        inla_results$spde,
        do.transform = TRUE
    )
    return(
        exp(rbind(
            result_field$summary.log.range.nominal,
            result_field$summary.log.variance.nominal
        ))
    )
}


extract_range_and_var <- function(params) {
    params <- params[c(5, 4, 6)]
    params$significant <- apply(params > 0, 1, function(x) var(x) == 0)
    params <- params %>%
        mutate(asterisk = case_when(
            significant ~ "*",
            !significant ~ ""
        ))
    params[1:3] <- params[1:3] %>% round(2)
    return(
        sprintf(
            "%s (%s)%s",
            params[[1]],
            paste(params[[2]], params[[3]], sep = "-"),
            params[[5]]
        )
    )
}


corr_plot <- function(df, var) {
    return(ggplot(data = df) +
        geom_point(aes(x = mean, y = error)) +
        theme_minimal() +
        xlab(var) +
        ylab("Average Forecast Error"))
}


inla_corr_plot <- function(ecoli_abs, staph_abs) {
    eco_spatial_effects <- lapply(ecoli_abs, load_spatial_field, "Ecoli")
    sau_spatial_effects <- lapply(staph_abs, load_spatial_field, "Staph")
    spatial_effects <- do.call(
        "rbind",
        c(eco_spatial_effects, sau_spatial_effects)
    )
    spatial_effects$param <- row.names(spatial_effects)
    spatial_effects <- spatial_effects[order(spatial_effects$param), ]
    spatial_effects$species <- c(rep("Ecoli", 4), rep("Staph", 4))
    spatial_effects$ab <- c(ecoli_abs, staph_abs)
    spatial_effects$param <- str_split(spatial_effects$param, "\\.", simplify = TRUE)[, 1]

    ecoli_results <- lapply(
        ecoli_abs,
        load_results,
        species = "Ecoli",
        inla_likelihood = "beta",
        percentage_correction = TRUE
    ) %>% add_ab_and_species(ecoli_abs, "Ecoli")
    staph_results <- lapply(
        staph_abs,
        load_results,
        species = "Staph",
        inla_likelihood = "beta",
        percentage_correction = TRUE
    ) %>% add_ab_and_species(staph_abs, "Staph")
    all_results <- rbind(ecoli_results, staph_results)
    results <- all_results %>%
        group_split(species, ab) %>%
        sapply(function(x) c(mean(x$error), x$species[[1]], x$ab[[1]])) %>%
        t() %>%
        as.data.frame()
    names(results) <- c("error", "species", "ab")

    df <- merge(results, spatial_effects)
    df$mean <- as.numeric(df$mean)
    df$error <- as.numeric(as.character(df$error))
    range_df <- df[df$param == "range", ]
    var_df <- df[df$param == "variance", ]

    p <- cowplot::plot_grid(
        corr_plot(range_df, "Range"),
        corr_plot(var_df, "Variance"),
        nrow = 1
    )
    ggsave("range_var_inla_accuracy_plot.png",
        plot = p
    )
}


ensemble_violinplot <- function(df) {
    best_other_forecast <- 
        order_forecast_accuracies(
            df[df$forecast_type != "ensemble_forecast",]
        )[[1]]
    df <- 
        df[
            df$forecast_type %in%
            c(best_other_forecast, "ensemble_forecast"),
        ]
    result <- format_results(df, c(best_other_forecast, "ensemble_forecast"))
    p <- ggplot(
        result,
        aes(
            x = forecast_type,
            y = error,
        )
    ) +
        geom_violin(draw_quantiles = 0.5) +
        theme_minimal() +
        theme(legend.position = "none") +
        xlab("Forecast Type") +
        ylab("Forecast Error")
    return(p)
}


ensemble_panel_plots <- function(ensemble_model = c("arima", "ets")) {
    get_results <- function(species, abs_list) {
        return(
            lapply(
                abs_list,
                function(x)
                    load_results(
                        species,
                        x,
                        inla_likelihood = "beta",
                        percentage_correction = TRUE,
                        ensemble_model = ensemble_model)
                    %>% mutate(ab = x)
            )
        )
    }
    results <- c(
        get_results("Ecoli", ecoli_abs),
        get_results("Staph", staph_abs)
    )
    plot_list <- lapply(results, ensemble_violinplot)
    panel_plot <- cowplot::plot_grid(
        plotlist = plot_list,
        ncol = 2,
        labels =  c(
            paste("Ecoli", ecoli_abs, sep = "-"),
            paste("Staph", staph_abs, sep = "-")
        ),
        label_x = 0.1,
        label_y = 1.1,
        label_size = 9
    ) +
        theme(plot.margin = margin(20))
    ggsave("forecast_ensemble_plot.png", plot = panel_plot)
}


sort_inla_accuracies <- function(species, ab) {
    inla_result <- load_inla_result(
        species,
        ab,
        likelihood = "beta",
        percentage_correction = TRUE
    )
    inla_result$error <- abs(
        inla_result$inla_forecast - inla_result$Perc.R
    )
    location_errors <- 
        inla_result %>%
            group_by(Location) %>%
            summarise_at(vars(error), mean)
    location_errors <- location_errors[order(location_errors$error), ]
    location_errors$ordinal_error <- as.numeric(as.factor(location_errors$error))
    location_errors$species <- species
    location_errors$ab <- ab
    return(location_errors)
}


plot_location_accuracies <- function(species, abs_list) {
    data <- do.call(
        "rbind",
        lapply(abs_list, function(x) sort_inla_accuracies(species, x))
    )
    wide_data <- pivot_wider(
        data,
        names_from = "ab",
        id_cols = "Location",
        values_from = "ordinal_error"
    )
    wide_data$Location <- as.factor(wide_data$Location)
    p <- ggparcoord(wide_data,
            columns = 2:5, groupColumn = 1, order = "anyClass",
            showPoints = TRUE,
            alphaLines = 0.3,
            title = paste(species, "error per location")
    ) +
    scale_color_viridis(discrete = TRUE) +
    xlab("Antibiotic") +
    ylab("Mean Absolute Error")
    ggsave(sprintf("%s_location_acc_plot.png", species), p)
}


panel_plots("Ecoli", ecoli_abs)
panel_plots("Staph", staph_abs)
tabulate_spatial_effect("Ecoli", ecoli_abs)
tabulate_spatial_effect("Staph", staph_abs)