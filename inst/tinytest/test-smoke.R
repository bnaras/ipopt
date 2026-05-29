v <- ipopt_version()
expect_equal(names(v), c("major", "minor", "release"))
expect_true(is.na(v[["major"]]) || v[["major"]] >= 3L)

ans <- ipopt_solve(
  x0 = 0,
  lower = -Inf,
  upper = Inf,
  eval_f = function(x) (x[1] - 2)^2,
  eval_grad_f = function(x) 2 * (x[1] - 2),
  hessian_structure = matrix(c(1L, 1L), ncol = 2),
  eval_h = function(x, obj_factor, lambda) 2 * obj_factor,
  options = list(print_level = 0L, tol = 1e-8)
)

expect_equal(ans$status_code, 0L)
expect_equal(ans$solution, 2, tolerance = 1e-6)
expect_equal(ans$objective, 0, tolerance = 1e-10)

ans <- ipopt_solve(
  x0 = c(0, 0),
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf),
  constraint_lower = 4,
  constraint_upper = 4,
  eval_f = function(x) sum((x - c(1, 2))^2),
  eval_grad_f = function(x) 2 * (x - c(1, 2)),
  eval_g = function(x) sum(x),
  jacobian_structure = matrix(c(1L, 1L, 1L, 2L), ncol = 2),
  eval_jac_g = function(x) c(1, 1),
  hessian_structure = matrix(c(1L, 1L, 2L, 2L), ncol = 2),
  eval_h = function(x, obj_factor, lambda) c(2 * obj_factor, 2 * obj_factor),
  options = list(print_level = 0L, tol = 1e-8)
)

expect_equal(ans$status_code, 0L)
expect_equal(ans$solution, c(1.5, 2.5), tolerance = 1e-6)
expect_equal(ans$constraints, 4, tolerance = 1e-8)
expect_equal(ans$objective, 0.5, tolerance = 1e-8)

# Linear objective with only bound constraints: the Lagrangian Hessian is
# structurally empty. In exact-Hessian mode (the default) this must still solve;
# the Hessian callback is supplied but reports zero entries.
ans <- ipopt_solve(
  x0 = c(0, 0),
  lower = c(1, 1),
  upper = c(2, 2),
  eval_f = function(x) sum(x),
  eval_grad_f = function(x) c(1, 1),
  hessian_structure = matrix(integer(), ncol = 2),
  eval_h = function(x, obj_factor, lambda) numeric(),
  options = list(print_level = 0L, tol = 1e-8)
)

expect_equal(ans$status_code, 0L)
expect_equal(ans$solution, c(1, 1), tolerance = 1e-6)
expect_equal(ans$objective, 2, tolerance = 1e-8)
