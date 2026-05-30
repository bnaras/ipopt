# ipopt

<!-- badges: start -->
  [![R-CMD-check](https://github.com/bnaras/ipopt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/ipopt/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

Lightweight R interface to the Ipopt nonlinear optimizer through Ipopt's C
interface.

The package links against the system Ipopt library, located at build time
with `pkg-config`. It is a `unix`-only package (Linux and macOS): there is
no Windows build, because Ipopt is not packaged for the Windows R toolchain.

## Installation

A system Ipopt (and `pkg-config`) must be available; the `configure` script
requires

```sh
pkg-config --cflags --libs ipopt
```

to succeed.

### Linux

Install the system library, then build from source:

```sh
# Debian / Ubuntu
sudo apt-get install -y coinor-libipopt-dev pkg-config
R CMD INSTALL .
```

A binary will also be available from the
[r-universe](https://bnaras.r-universe.dev) once Ipopt is registered in the
system-requirements database:

```r
install.packages("ipopt", repos = "https://bnaras.r-universe.dev")
```

### macOS (build from source)

Homebrew provides Ipopt, so macOS users build from source. The `configure`
script locates the Homebrew Ipopt automatically:

```sh
brew install ipopt pkg-config
```

```r
install.packages("ipopt", repos = "https://bnaras.r-universe.dev", type = "source")
```

Equivalently, clone the repository and run `R CMD INSTALL .`.

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
