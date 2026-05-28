#' Ipopt Library Version
#'
#' @return An integer vector with elements `major`, `minor`, and `release`.
#' @export
ipopt_version <- function() {
  .Call("C_ipopt_version")
}

#' Solve a Nonlinear Program with Ipopt
#'
#' This is a low-level wrapper over Ipopt's C interface. Derivative callbacks
#' are supplied as R functions. Indices in `jacobian_structure` and
#' `hessian_structure` are one-based R indices and are converted to Ipopt's
#' zero-based C index style internally.
#'
#' @param x0 Numeric starting point.
#' @param lower,upper Numeric variable bounds, length `length(x0)`.
#' @param constraint_lower,constraint_upper Numeric constraint bounds.
#' @param eval_f Function `f(x)` returning the scalar objective value.
#' @param eval_grad_f Function `grad_f(x)` returning a numeric vector.
#' @param eval_g Function `g(x)` returning a numeric vector of constraints.
#' @param eval_jac_g Function returning the Jacobian nonzero values in the
#'   order specified by `jacobian_structure`.
#' @param jacobian_structure Integer two-column matrix of one-based row/column
#'   positions for Jacobian nonzeros.
#' @param eval_h Function `(x, obj_factor, lambda)` returning Hessian nonzero
#'   values in the order specified by `hessian_structure`.
#' @param hessian_structure Integer two-column matrix of one-based lower
#'   triangular row/column positions for Hessian nonzeros.
#' @param options Named list of Ipopt options. Numeric, integer, logical, and
#'   character values are supported.
#'
#' @return A list with solution, objective, constraints, multipliers, and status.
#' @export
ipopt_solve <- function(x0,
                        lower,
                        upper,
                        constraint_lower = numeric(),
                        constraint_upper = numeric(),
                        eval_f,
                        eval_grad_f,
                        eval_g = NULL,
                        eval_jac_g = NULL,
                        jacobian_structure = matrix(integer(), ncol = 2),
                        eval_h = NULL,
                        hessian_structure = matrix(integer(), ncol = 2),
                        options = list()) {
  stopifnot(is.numeric(x0), is.numeric(lower), is.numeric(upper))
  n <- length(x0)
  if (length(lower) != n || length(upper) != n) {
    stop("`lower` and `upper` must have the same length as `x0`.", call. = FALSE)
  }
  m <- length(constraint_lower)
  if (length(constraint_upper) != m) {
    stop("constraint bounds must have the same length.", call. = FALSE)
  }
  if (m > 0L && (is.null(eval_g) || is.null(eval_jac_g))) {
    stop("constraint callbacks are required when constraints are present.", call. = FALSE)
  }
  jacobian_structure <- .as_structure(jacobian_structure, "jacobian_structure")
  hessian_structure <- .as_structure(hessian_structure, "hessian_structure")

  .Call(
    "C_ipopt_solve",
    as.numeric(x0),
    as.numeric(lower),
    as.numeric(upper),
    as.numeric(constraint_lower),
    as.numeric(constraint_upper),
    eval_f,
    eval_grad_f,
    eval_g %||% function(x) numeric(),
    eval_jac_g %||% function(x) numeric(),
    jacobian_structure,
    eval_h %||% function(x, obj_factor, lambda) numeric(),
    hessian_structure,
    options
  )
}

.as_structure <- function(x, name) {
  x <- as.matrix(x)
  if (!is.integer(x)) storage.mode(x) <- "integer"
  if (!identical(ncol(x), 2L)) {
    stop(sprintf("`%s` must have two columns.", name), call. = FALSE)
  }
  x
}

`%||%` <- function(x, y) if (is.null(x)) y else x
