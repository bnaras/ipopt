# ipopt

Lightweight R interface to the Ipopt nonlinear optimizer through Ipopt's C
interface.

This package is intended first for r-universe builds where Ipopt and its open
linear solver dependencies, such as MUMPS, are available as system libraries
discoverable with `pkg-config`.

## Build

```sh
R CMD INSTALL .
```

The package configure script requires:

```sh
pkg-config --cflags --libs ipopt
```

to succeed.

## Smoke Test

```r
library(ipopt)

ipopt_solve(
  x0 = 0,
  lower = -Inf,
  upper = Inf,
  eval_f = function(x) (x[1] - 2)^2,
  eval_grad_f = function(x) 2 * (x[1] - 2),
  hessian_structure = matrix(c(1L, 1L), ncol = 2),
  eval_h = function(x, obj_factor, lambda) 2 * obj_factor,
  options = list(print_level = 0L)
)
```
