rosen <- function(x, y) {
  (1 - x)^2 + 100 * (y - x^2)^2
}

rosen.gr <- function(x, y) {
  c(-2 * (1 - x) - 400 * x * (y - x^2), 200 * (y - x^2))
}

rosen.bfgs <- function(par, control = list()) {
  control <- modifyList(list(maxit = 100, gtol = 1e-6), control)
  
  bfgs(par,
       function(w){rosen(w[1], w[2])},
       function(w){rosen.gr(w[1], w[2])}, control = control)
}

result <- rosen.bfgs(c(-1.2, 1), control = list(maxit = 100, gtol = 1e-7))
