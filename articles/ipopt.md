# Solving nonlinear programs with ipopt

## The problem Ipopt solves

[Ipopt](https://github.com/coin-or/Ipopt) (Interior Point OPTimizer) is
an open-source solver for large-scale **nonlinear programs** of the form

``` math
\begin{aligned}
\min_{x \in \mathbb{R}^n} \quad & f(x) \\
\text{subject to} \quad & g^{\mathrm{L}} \le g(x) \le g^{\mathrm{U}}, \\
                        & x^{\mathrm{L}} \le x \le x^{\mathrm{U}},
\end{aligned}
```

where $`f : \mathbb{R}^n \to \mathbb{R}`$ is the objective,
$`g : \mathbb{R}^n \to
\mathbb{R}^m`$ collects the (possibly nonlinear) constraints, and the
bounds may be infinite. Equality constraints are expressed by setting
the corresponding lower and upper bounds equal; one-sided inequalities
use an infinite bound on the open side.

This package is a thin wrapper over Ipopt’s C interface. You describe
the problem by supplying R functions that evaluate $`f`$, its gradient,
the constraints, and the relevant derivative matrices. The single entry
point is
[`ipopt_solve()`](https://bnaras.github.io/ipopt/reference/ipopt_solve.md).

## Anatomy of `ipopt_solve()`

``` r

ipopt_solve(
  x0,                                       # starting point, length n
  lower, upper,                             # variable bounds x^L, x^U (length n)
  constraint_lower = numeric(),             # constraint bounds g^L, g^U (length m)
  constraint_upper = numeric(),
  eval_f, eval_grad_f,                      # objective + its gradient
  eval_g = NULL, eval_jac_g = NULL,         # constraints + their Jacobian values
  jacobian_structure = matrix(integer(), ncol = 2),
  eval_h = NULL,                            # Hessian of the Lagrangian (values)
  hessian_structure = matrix(integer(), ncol = 2),
  options = list()                          # named list of Ipopt options
)
```

The interesting part — and what this vignette is really about — is how
the **callbacks** and their **sparsity structures** fit together.

### Bounds and constraints

- `lower` / `upper` are the variable bounds, each of length `n`. Use
  `-Inf` / `Inf` for a free variable, and set `lower[i] == upper[i]` to
  fix variable `i`.

- `constraint_lower` / `constraint_upper` are the bounds on $`g(x)`$,
  each of length `m`. The flavour of each constraint is encoded entirely
  in its bounds:

  | Constraint             | `constraint_lower` | `constraint_upper` |
  |------------------------|:------------------:|:------------------:|
  | $`g_j(x) = b`$         |        `b`         |        `b`         |
  | $`g_j(x) \ge a`$       |        `a`         |       `Inf`        |
  | $`g_j(x) \le c`$       |       `-Inf`       |        `c`         |
  | $`a \le g_j(x) \le c`$ |        `a`         |        `c`         |

When there are no constraints (`m == 0`) you omit `eval_g`,
`eval_jac_g`, and `jacobian_structure` entirely.

### Objective: `eval_f` and `eval_grad_f`

``` r

eval_f      <- function(x) ...      # returns a scalar, f(x)
eval_grad_f <- function(x) ...      # returns the dense length-n gradient ∇f(x)
```

The gradient is **dense**: return all `n` partial derivatives in
variable order, even the zeros.

### Constraints: `eval_g`, `eval_jac_g`, and `jacobian_structure`

``` r

eval_g     <- function(x) ...       # returns the length-m vector g(x)
eval_jac_g <- function(x) ...       # returns the *nonzero* Jacobian values
```

The constraint Jacobian $`J_{jk} = \partial g_j / \partial x_k`$ is an
$`m \times n`$ matrix that is usually sparse, so Ipopt asks for it in
coordinate form. You declare the **sparsity pattern once** in
`jacobian_structure`, a two-column integer matrix whose rows are the
`(row, col)` positions of the nonzeros, using **one-based** R indices:

``` r

# row j, column k of the Jacobian:
jacobian_structure <- cbind(row = c(...), col = c(...))
```

`eval_jac_g(x)` then returns just the nonzero **values**, in the *same
row order* as `jacobian_structure`. Ipopt evaluates the values many
times during the solve, but the structure is fixed and read only once.
(Internally the one-based indices are translated to Ipopt’s zero-based C
convention for you.)

### Hessian: `eval_h` and `hessian_structure`

This is the step most people get wrong. Ipopt does not want the Hessian
of the objective — it wants the Hessian of the **Lagrangian**

``` math
\nabla^2_{xx} L(x, \sigma, \lambda)
  = \sigma\, \nabla^2 f(x) \;+\; \sum_{j=1}^{m} \lambda_j\, \nabla^2 g_j(x),
```

a symmetric $`n \times n`$ matrix. Because it is symmetric, you supply
**only its lower triangle** (entries with `row >= col`):

``` r

# lower-triangular nonzeros, one-based, row >= col:
hessian_structure <- cbind(row = c(...), col = c(...))

eval_h <- function(x, obj_factor, lambda) ...
```

- `obj_factor` is the scalar $`\sigma`$ multiplying the objective
  Hessian.
- `lambda` is the length-`m` vector of constraint multipliers
  $`\lambda`$.
- The return value is the nonzero lower-triangular values, in the same
  row order as `hessian_structure`.

If your constraints are linear their Hessians vanish and only
$`\sigma \nabla^2 f(x)`$ remains; if the objective is linear too, the
Lagrangian Hessian is structurally empty (a valid case — pass an empty
`hessian_structure` and an `eval_h` that returns
[`numeric()`](https://rdrr.io/r/base/numeric.html)).

**Skipping the exact Hessian.** If a second derivative is inconvenient,
pass `eval_h = NULL` together with a quasi-Newton approximation:

``` r

options = list(hessian_approximation = "limited-memory")
```

Ipopt then builds an L-BFGS approximation from gradients alone. (Passing
`eval_h = NULL` *without* requesting an approximation is an error, since
the exact mode has nothing to evaluate.)

### Options

`options` is a named list passed straight to Ipopt. Values that are
character, integer/logical, or numeric are routed to the string,
integer, and number option setters respectively. A few useful ones:

| Option                  | Meaning                                   |
|-------------------------|-------------------------------------------|
| `print_level`           | 0 (silent) … 12 (verbose); default 5      |
| `sb = "yes"`            | suppress the startup banner               |
| `tol`                   | convergence tolerance                     |
| `max_iter`              | iteration cap                             |
| `hessian_approximation` | `"exact"` (default) or `"limited-memory"` |
| `mu_strategy`           | `"monotone"` or `"adaptive"`              |

The full list is in the [Ipopt options
reference](https://coin-or.github.io/Ipopt/OPTIONS.html).

### The return value

[`ipopt_solve()`](https://bnaras.github.io/ipopt/reference/ipopt_solve.md)
returns a list:

| Field                     | Contents                             |
|---------------------------|--------------------------------------|
| `solution`                | the optimal $`x^\star`$              |
| `objective`               | $`f(x^\star)`$                       |
| `constraints`             | $`g(x^\star)`$                       |
| `status_code`             | Ipopt return code (`0` = success)    |
| `constraint_multipliers`  | $`\lambda`$ at the solution          |
| `lower_bound_multipliers` | bound multipliers $`z^{\mathrm{L}}`$ |
| `upper_bound_multipliers` | bound multipliers $`z^{\mathrm{U}}`$ |

A `status_code` of `0` means *Solve_Succeeded*; `1` is
*Solved_To_Acceptable_Level*. Negative codes signal trouble (e.g. `-1`
maximum iterations, `-2` restoration failed). The constraint multipliers
are the Lagrange multipliers $`\lambda`$ — at the solution they carry
the usual sensitivity (shadow-price) interpretation.

We silence the solver in the examples below with:

``` r

quiet <- list(print_level = 0L, sb = "yes")
```

## Example 1: HS071 — the full machinery

The Hock–Schittkowski problem \#71 is the canonical Ipopt tutorial
example because it exercises everything at once: a nonlinear objective,
one nonlinear inequality, one nonlinear equality, and bounds.

``` math
\begin{aligned}
\min_{x \in \mathbb{R}^4} \quad & x_1 x_4 (x_1 + x_2 + x_3) + x_3 \\
\text{s.t.}\quad & x_1 x_2 x_3 x_4 \ge 25, \\
                 & x_1^2 + x_2^2 + x_3^2 + x_4^2 = 40, \\
                 & 1 \le x_i \le 5, \quad i = 1,\dots,4.
\end{aligned}
```

**Gradient.** With $`f = x_1 x_4(x_1 + x_2 + x_3) + x_3`$,

``` math
\nabla f = \begin{pmatrix}
x_4(2x_1 + x_2 + x_3) \\ x_1 x_4 \\ x_1 x_4 + 1 \\ x_1(x_1 + x_2 + x_3)
\end{pmatrix}.
```

**Jacobian.** Both constraints depend on every variable, so the Jacobian
is dense ($`2 \times 4 = 8`$ nonzeros). Row 1 is
$`\nabla(x_1 x_2 x_3 x_4)`$, row 2 is $`\nabla(\sum x_i^2)`$.

**Lagrangian Hessian.**
$`\nabla^2 L = \sigma \nabla^2 f + \lambda_1 \nabla^2 g_1
+ \lambda_2 \nabla^2 g_2`$ is a full symmetric $`4 \times 4`$ matrix;
its lower triangle has 10 nonzeros. We assemble the lower triangle row
by row.

``` r

hs071 <- ipopt_solve(
  x0    = c(1, 5, 5, 1),
  lower = rep(1, 4),
  upper = rep(5, 4),
  constraint_lower = c(25, 40),     # g1 >= 25  (upper Inf); g2 == 40
  constraint_upper = c(Inf, 40),

  eval_f = function(x) x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3],

  eval_grad_f = function(x) c(
    x[4] * (2 * x[1] + x[2] + x[3]),
    x[1] * x[4],
    x[1] * x[4] + 1,
    x[1] * (x[1] + x[2] + x[3])
  ),

  eval_g = function(x) c(prod(x), sum(x^2)),

  # dense 2 x 4 Jacobian: rows (1,1,1,1, 2,2,2,2), cols (1,2,3,4, 1,2,3,4)
  jacobian_structure = cbind(rep(1:2, each = 4), rep(1:4, 2)),
  eval_jac_g = function(x) c(
    x[2] * x[3] * x[4], x[1] * x[3] * x[4],          # d g1 / d x1, x2
    x[1] * x[2] * x[4], x[1] * x[2] * x[3],          # d g1 / d x3, x4
    2 * x[1], 2 * x[2], 2 * x[3], 2 * x[4]           # d g2 / d x
  ),

  # lower triangle of the 4 x 4 Lagrangian Hessian (row >= col), row by row
  hessian_structure = cbind(
    c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4),
    c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4)
  ),
  eval_h = function(x, obj_factor, lambda) {
    of <- obj_factor; l1 <- lambda[1]; l2 <- lambda[2]
    c(
      of * (2 * x[4])                  + l2 * 2,        # (1,1)
      of * x[4]            + l1 * (x[3] * x[4]),        # (2,1)
      0                                + l2 * 2,        # (2,2)
      of * x[4]            + l1 * (x[2] * x[4]),        # (3,1)
                             l1 * (x[1] * x[4]),        # (3,2)
      0                                + l2 * 2,        # (3,3)
      of * (2 * x[1] + x[2] + x[3]) + l1 * (x[2] * x[3]),  # (4,1)
      of * x[1]            + l1 * (x[1] * x[3]),        # (4,2)
      of * x[1]            + l1 * (x[1] * x[2]),        # (4,3)
      0                                + l2 * 2         # (4,4)
    )
  },

  options = c(quiet, tol = 1e-8)
)

hs071$status_code
#> [1] 0
hs071$solution
#> [1] 1.000000 4.743000 3.821150 1.379408
hs071$objective
#> [1] 17.01402
hs071$constraints                 # should hit 25 and 40
#> [1] 25 40
```

The optimum is $`x^\star \approx (1, 4.743, 3.821, 1.379)`$ with
objective $`\approx 17.014`$, the textbook answer.

## Example 2: Rosenbrock — unconstrained, exact Hessian

The Rosenbrock “banana” function is the classic unconstrained test
problem:

``` math
\min_{x, y}\; (1 - x)^2 + 100\,(y - x^2)^2,
```

with the global minimum $`0`$ at $`(1, 1)`$ sitting at the bottom of a
curved valley. With no constraints there is no Jacobian, and the
Lagrangian Hessian is just $`\sigma \nabla^2 f`$ (so `lambda` is
unused). The Hessian’s lower triangle has three entries: $`(1,1)`$,
$`(2,1)`$, $`(2,2)`$.

``` r

rosen <- ipopt_solve(
  x0    = c(-1.2, 1),               # the traditional hard starting point
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf),

  eval_f = function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2,

  eval_grad_f = function(x) c(
    -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2)
  ),

  hessian_structure = cbind(c(1, 2, 2), c(1, 1, 2)),
  eval_h = function(x, obj_factor, lambda) obj_factor * c(
    2 - 400 * x[2] + 1200 * x[1]^2,   # (1,1)
    -400 * x[1],                      # (2,1)
    200                               # (2,2)
  ),

  options = c(quiet, tol = 1e-10)
)

rosen$status_code
#> [1] 0
rosen$solution
#> [1] 1 1
rosen$objective
#> [1] 0
```

## Example 3: Rosenbrock without a Hessian (limited-memory)

Coding (and debugging) an exact Lagrangian Hessian is the most
error-prone part of using Ipopt. When you would rather not, drop
`eval_h` and let Ipopt build an L-BFGS approximation from gradients. The
objective and gradient are identical to Example 2; everything
Hessian-related disappears:

``` r

rosen_lm <- ipopt_solve(
  x0    = c(-1.2, 1),
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf),

  eval_f = function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2,
  eval_grad_f = function(x) c(
    -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2)
  ),

  eval_h = NULL,
  options = c(quiet, tol = 1e-8, hessian_approximation = "limited-memory")
)

rosen_lm$status_code
#> [1] 0
rosen_lm$solution
#> [1] 1 1
```

Same optimum, no second derivatives. The cost is typically more
iterations; the benefit is far less code to get wrong.

## Example 4: Maximum-entropy distribution — equality constraints

A staple of statistics and information-theory courses: among all
probability distributions $`p`$ on $`\{1, \dots, K\}`$ with a prescribed
mean $`\mu`$, find the one of **maximum entropy**. Equivalently we
minimize the negative entropy

``` math
\min_{p \in \mathbb{R}^K}\; \sum_{i=1}^{K} p_i \log p_i
\quad\text{s.t.}\quad
\sum_i p_i = 1, \qquad \sum_i i\,p_i = \mu, \qquad p_i \ge 0.
```

Both constraints are **linear**, so their Hessians vanish and the
Lagrangian Hessian is just
$`\sigma \nabla^2 f = \sigma\,\mathrm{diag}(1/p_i)`$ — a **diagonal**
matrix, so `hessian_structure` lists only the diagonal. The gradient is
$`\partial f/\partial p_i = \log p_i + 1`$, and the (constant) Jacobian
rows are $`\mathbf{1}`$ and $`(1, 2, \dots, K)`$.

``` r

K       <- 5
support <- 1:K
mu      <- 2                       # desired mean (uniform mean would be 3)

ent <- ipopt_solve(
  x0    = rep(1 / K, K),            # start at the uniform distribution
  lower = rep(1e-8, K),            # keep p_i > 0 so log p_i is finite
  upper = rep(1, K),
  constraint_lower = c(1, mu),     # both equalities: lower == upper
  constraint_upper = c(1, mu),

  eval_f      = function(p) sum(p * log(p)),
  eval_grad_f = function(p) log(p) + 1,

  eval_g = function(p) c(sum(p), sum(support * p)),
  jacobian_structure = cbind(rep(1:2, each = K), rep(1:K, 2)),
  eval_jac_g = function(p) c(rep(1, K), support),

  hessian_structure = cbind(1:K, 1:K),               # diagonal only
  eval_h = function(p, obj_factor, lambda) obj_factor * (1 / p),

  options = c(quiet, tol = 1e-9)
)

ent$status_code
#> [1] 0
round(ent$solution, 5)
#> [1] 0.45936 0.26079 0.14806 0.08406 0.04772
c(total = sum(ent$solution), mean = sum(support * ent$solution))
#> total  mean 
#>     1     2
```

Maximum-entropy theory says the answer must be a discrete exponential
(Gibbs) distribution $`p_i \propto e^{-\theta i}`$, and the multiplier
on the mean constraint *is* that natural parameter $`\theta`$. We can
confirm the solver recovered exactly that:

``` r

theta <- ent$constraint_multipliers[2]
closed_form <- exp(-theta * support)
closed_form <- closed_form / sum(closed_form)
round(closed_form, 5)              # matches ent$solution
#> [1] 0.45936 0.26079 0.14806 0.08406 0.04772
```

The numerical solution and the closed form agree to solver tolerance — a
nice illustration that the `constraint_multipliers` returned by
[`ipopt_solve()`](https://bnaras.github.io/ipopt/reference/ipopt_solve.md)
are the genuine Lagrange multipliers, not a by-product.

## Where to go next

- [`?ipopt_solve`](https://bnaras.github.io/ipopt/reference/ipopt_solve.md)
  for the argument-level reference.
- The [Ipopt documentation](https://coin-or.github.io/Ipopt/) for the
  full options catalogue and the theory behind the interior-point
  method.
- For modelling languages that *generate* these callbacks automatically,
  see the [CVXR](https://cvxr.rbind.io) project, which can target Ipopt
  as a backend for smooth nonlinear programs.
