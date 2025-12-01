#' @title transform X
#'
#' @description Take data vector X and convert factors with numeric. Binary
#'   factors will be encoded by {0, 1}. Factors with more than two levels will
#'   be transformed using one-hot encoding, dropping one level to avoud
#'   multicollinearity.
#'
#' @param x Design matrix X with covariates
#'
#' @return Transformed design matrix X
#'
#' @author Peter Degen
#'
#' @export
transform_X <- function(X) {
  X <- as.data.frame(X)
  out_list <- list()

  for (nm in names(X)) {
    col <- X[[nm]]

    if (is.factor(col)) {
      lev <- levels(col)
      k <- length(lev)

      if (k == 2) {
        # Binary factor -> convert to 0/1 and keep original column name
        out_list[[nm]] <- as.numeric(col) - 1

      } else {
        # Multi-level factor -> one-hot encode, drop first level
        mm <- model.matrix(~ col)[, -1, drop = FALSE]

        # Rename dummy columns as col_level
        new_names <- paste0(nm, "_", lev[-1])
        colnames(mm) <- new_names

        out_list[[nm]] <- mm
      }

    } else {
      # Non-factor: numeric or character (character would error)
      out_list[[nm]] <- col
    }
  }

  # Combine results
  X_mat <- do.call(cbind, out_list)
  X_mat <- as.matrix(X_mat)
  storage.mode(X_mat) <- "numeric"
  colnames(X_mat) <- names(X)
  X_mat
}