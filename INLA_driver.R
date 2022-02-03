library(INLA)

build_mesh_and_spde <- function(data, x_col, y_col, time_col) {
    locations_matrix <- cbind(
        unique(data[[x_col]]),
        unique(data[[y_col]])
    )
    mesh <- INLA::inla.mesh.2d(
        loc = locations_matrix, max.edge = 0.2, offset = 0.2, min.angle = 1
    )
    return(inla.spde2.matern(mesh = mesh))
}


build_stack <- function(data, validation_data, x_col,
                        y_col, time_col, response_col, spde) {
    n_pred_days <- length(unique(data[[time_col]]))
    n_val_days <- length(unique(validation_data[[time_col]]))

    a_est <-
        inla.spde.make.A(spde,
            loc = as.matrix(data[c(x_col, y_col)]),
            group = data[[time_col]],
            n.group = n_pred_days
        )
    est_field_indices <-
        inla.spde.make.index("field",
            n.spde = spde$n.spde,
            n.group = n_pred_days
        )
    stack_est <- inla.stack(
        data = list(outcome = data[[response_col]]),
        A = list(a_est, 1),
        effects = list(
            c(
                est_field_indices,
                list(Intercept = 1)
            ),
            list(data[c(x_col, y_col, time_col)])
        ),
        tag = "est"
    )

    a_pred <-
        inla.spde.make.A(spde,
            loc = as.matrix(validation_data[c(x_col, y_col)]),
            group = as.numeric(as.factor(validation_data[[time_col]])),
            n.group = n_val_days
        )
    pred_field_indices <-
        inla.spde.make.index("field",
            n.spde = spde$n.spde,
            n.group = n_val_days
        )
    stack_pred <- inla.stack(
        data = list(outcome = NA),
        A = list(a_pred, 1),
        effects = list(
            c(
                pred_field_indices,
                list(Intercept = 1)
            ),
            list(validation_data[c(x_col, y_col, time_col)])
        ),
        tag = "pred"
    )

    return(list(stack_est, stack_pred))
}