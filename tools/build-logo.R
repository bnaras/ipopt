## Generate the ipopt hex sticker.
##
## Source of truth = THIS script. The wireframe is the Rosenbrock "banana"
## valley (z = (1 - x)^2 + 100 (y - x^2)^2), the canonical nonlinear test
## function Ipopt is built to minimize, generated in R so colour, view and
## density are fully under our control and the mesh stays contained inside the
## hex. Renders man/figures/logo.png (+ logo@2x.png for Retina) and refreshes
## the pkgdown favicons. Run from the package root:
##
##   Rscript tools/build-logo.R
##
## This script is .Rbuildignore'd, so it is not shipped in the tarball
## (tools/version.c, used by configure, is kept).
##
## Palette (COIN-OR teal family):
##   Hex fill ...... #115e59  deep teal
##   Hex border .... #0b3b38  darker teal
##   Mesh .......... #c7efe9  pale aqua wireframe
##   Wordmark ...... #ffffff  white "ipopt"

suppressMessages({library(hexSticker); library(grDevices)})

TEAL   <- "#115e59"
BORDER <- "#0b3b38"
MESH   <- "#c7efe9"
WORD   <- "#ffffff"

## Rosenbrock wireframe -> transparent PNG. The surface is clamped (sqrt-scaled
## height) so the steep walls do not blow up the persp box; what reads is the
## characteristic curved valley sweeping across the cell.
mesh_png <- function(file, col = MESH, lwd = 1.0, theta = -32, phi = 26, n = 30) {
  x <- seq(-1.5, 1.5, length.out = n)
  y <- seq(-0.6, 1.7, length.out = n)
  z <- outer(x, y, function(x, y) (1 - x)^2 + 100 * (y - x^2)^2)
  z <- sqrt(z)                         # tame the dynamic range for the view
  grDevices::png(file, width = 1100, height = 1100, bg = "transparent")
  op <- par(mar = c(0, 0, 0, 0))
  persp(x, y, z, theta = theta, phi = phi, col = NA, border = col, lwd = lwd,
        box = FALSE, axes = FALSE, scale = TRUE, expand = 0.5, r = 3, d = 2)
  par(op); grDevices::dev.off(); file
}

build <- function(out, dpi) {
  m <- mesh_png(tempfile(fileext = ".png"))
  sticker(
    # valley mesh in the lower-centre; lower-case wordmark up top
    subplot = m, s_x = 1.0, s_y = 0.78, s_width = 0.70, s_height = 0.70,
    package = "ipopt", p_x = 1.0, p_y = 1.44, p_size = 30, p_color = WORD,
    p_family = "sans", p_fontface = "bold",
    h_fill = TEAL, h_color = BORDER, h_size = 1.4,
    dpi = dpi, filename = out
  )
  cat("wrote", out, "\n")
}

if (!dir.exists("man/figures")) dir.create("man/figures", recursive = TRUE)
build("man/figures/logo.png",    320)   # standard
build("man/figures/logo@2x.png", 640)   # Retina

if (requireNamespace("pkgdown", quietly = TRUE)) {
  cat("Rebuilding pkgdown favicons...\n")
  try(pkgdown::build_favicons(overwrite = TRUE), silent = TRUE)
}
