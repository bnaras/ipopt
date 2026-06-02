# ipopt

Lightweight R interface to the Ipopt nonlinear optimizer through Ipopt’s
C interface.

See the [*Solving nonlinear programs with
ipopt*](https://bnaras.github.io/ipopt/articles/ipopt.html) vignette for
a guided tour of the objective, gradient, constraint Jacobian, and
Hessian-of-the-Lagrangian callbacks, with worked examples (HS071,
Rosenbrock, limited-memory, and a maximum-entropy problem).

The package links against an Ipopt library located at build time. On
Linux and macOS this is the system Ipopt found with `pkg-config`; on
Windows it is one of Ipopt’s official prebuilt release binaries (see
below). There is no binary on CRAN/r-universe for every platform yet, so
most users build from source against a locally installed Ipopt.

## Installation

A system Ipopt (and `pkg-config`) must be available; the `configure`
script requires

``` sh
pkg-config --cflags --libs ipopt
```

to succeed.

### Linux

Install the system library, then build from source:

``` sh
# Debian / Ubuntu
sudo apt-get install -y coinor-libipopt-dev pkg-config
R CMD INSTALL .
```

A binary will also be available from the
[r-universe](https://bnaras.r-universe.dev) once Ipopt is registered in
the system-requirements database:

``` r

install.packages("ipopt", repos = "https://bnaras.r-universe.dev")
```

### macOS (build from source)

Homebrew provides Ipopt, so macOS users build from source. The
`configure` script locates the Homebrew Ipopt automatically:

``` sh
brew install ipopt pkg-config
```

``` r

install.packages("ipopt", repos = "https://bnaras.r-universe.dev", type = "source")
```

Equivalently, clone the repository and run `R CMD INSTALL .`.

### Windows (prebuilt Ipopt binary)

Ipopt does not ship in the Windows R toolchain, but the Ipopt project
publishes prebuilt Windows binaries that this package can link against.
Although those binaries are built with MSVC, they export the Ipopt C
interface as plain C symbols, so Rtools’ (MinGW) linker can use them
directly.

1.  Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
    (matching your R version).

2.  Download the latest **`win64-msvs2022-md`** zip from [Ipopt
    releases](https://github.com/coin-or/Ipopt/releases) and extract it,
    for example to `C:\Ipopt`, giving you `C:\Ipopt\include\coin-or`,
    `C:\Ipopt\lib`, and `C:\Ipopt\bin`.

3.  Build from source, pointing the package at that folder (use
    **forward slashes**):

    ``` sh
    R CMD INSTALL ipopt ^
      --configure-vars="INCLUDE_DIR=C:/Ipopt/include/coin-or LIB_DIR=C:/Ipopt/lib"
    ```

    Equivalently, set `IPOPT_HOME=C:/Ipopt` in the environment before
    installing and the script will find the headers, library, and DLL
    itself.

4.  At **run time**, the package loads `ipopt-3.dll` and its
    dependencies (`coinmumps-3.dll` and the Intel runtime DLLs, all in
    `C:\Ipopt\bin`). Make sure that folder is on the `PATH` for the R
    session, e.g.

    ``` r

    Sys.setenv(PATH = paste("C:/Ipopt/bin", Sys.getenv("PATH"), sep = ";"))
    library(ipopt)
    ```

    (Add `C:\Ipopt\bin` to the system `PATH` to avoid doing this each
    session.)

## Smoke Test

``` r

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
