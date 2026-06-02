# Solve a Nonlinear Program with Ipopt

This is a low-level wrapper over Ipopt's C interface. Derivative
callbacks are supplied as R functions. Indices in `jacobian_structure`
and `hessian_structure` are one-based R indices and are converted to
Ipopt's zero-based C index style internally.

## Usage

``` r
ipopt_solve(
  x0,
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
  options = list()
)
```

## Arguments

- x0:

  Numeric starting point.

- lower, upper:

  Numeric variable bounds, length `length(x0)`.

- constraint_lower, constraint_upper:

  Numeric constraint bounds.

- eval_f:

  Function `f(x)` returning the scalar objective value.

- eval_grad_f:

  Function `grad_f(x)` returning a numeric vector.

- eval_g:

  Function `g(x)` returning a numeric vector of constraints.

- eval_jac_g:

  Function returning the Jacobian nonzero values in the order specified
  by `jacobian_structure`.

- jacobian_structure:

  Integer two-column matrix of one-based row/column positions for
  Jacobian nonzeros.

- eval_h:

  Function `(x, obj_factor, lambda)` returning Hessian nonzero values in
  the order specified by `hessian_structure`. A structurally empty
  Hessian (no rows in `hessian_structure`, the function returning
  [`numeric()`](https://rdrr.io/r/base/numeric.html)) is valid and
  corresponds to a linear objective with linear/bound constraints. Pass
  `NULL` only when no exact Hessian is available, in which case a
  Hessian approximation must be requested via `options` (for example
  `hessian_approximation = "limited-memory"`); `eval_h = NULL` without
  such an option is an error.

- hessian_structure:

  Integer two-column matrix of one-based lower triangular row/column
  positions for Hessian nonzeros.

- options:

  Named list of Ipopt options. Numeric, integer, logical, and character
  values are supported.

## Value

A list with solution, objective, constraints, multipliers, and status.
