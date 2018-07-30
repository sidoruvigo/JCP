#' @title comp2regr.onesided
#' @description function that calculates the p-value based on c-v bandwidths
#' @export
comp2regr.onesided <- function(x1, x2, y1, y2, p = function(x){0.5}, w1 = function(x){1.}, w2 = function(x){1.}){
  # GENERAL FUNCTIONS (kernel-based estimators)

  # Kernel Epanechnikov (density)
  kernel <- function(u) {
    (0.75 * (1 - u ^ 2)) * (u < 1) * (u > -1)
  }

  # Kernel density estimator. This function calculates the kernel density estimator on the collection of "npoints" points "points" based on the "ndata" observations xdata and the bandwidth "h".
  kerneldensity <- function(ndata, data, points, h) {
    rowSums(kernel(outer(points, data, "-") / h)) / (ndata * h)
  }

  # Nadaraya-Watson estimator. This function calculates the N-W estimator on the collection of "npoints" points "points" based on the "ndata" observations (xdata,ydata) and the bandwidth "h".
  nadarayawatson <- function(ndata, xdata, ydata, npoints, points, h) {
    as.vector({matk <- kernel((points %*% t(rep(1, ndata)) - t(xdata %*% t(rep(1, npoints)))) / h)
    (matk %*% ydata) / (matk %*% rep(1, ndata))
    })
  }

  # Cross-validation bandwidth selector (Nadaraya-Watson)
  h.crossvalidation.nw <- function(n, x, y, hmin, hmax, hngrid) {
    crossvalue <- rep(0, hngrid)
    h <- seq(hmin, hmax, len = hngrid)
    for (j in 1:hngrid) {
      for (i in 1:n) {
        crossvalue[j] <- crossvalue[j] + (y[i] - nadarayawatson(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^
          2
      }
    }
    crossvalue <- crossvalue / n
    #plot(h,crossvalue)
    h.crossvalidation.nw <- h[which.min(crossvalue)]
  }

  # Local-linear regression estimator
  locallinear <- function(ndata, xdata, ydata, npoints, points, h) {
    # See Wand and Jones (1995). "Kernel Smoothing" (page 119).
    mat.x <- outer(points, xdata, "-")
    mat.k <- kernel(mat.x / h)
    mat.y <- matrix(rep(ydata, npoints), nrow = npoints, byrow = T)
    s0 <- matrix(rep(rowSums(mat.k) / ndata, ndata), ncol = ndata)
    s1 <- matrix(rep(rowSums(mat.x * mat.k) / ndata, ndata), ncol = ndata)
    s2 <- matrix(rep(rowSums((mat.x ^ 2) * mat.k) / ndata, ndata), ncol = ndata)
    return((rowSums((s2 - s1 * mat.x) * mat.k * mat.y) / ndata) / (s2[, 1] * s0[, 1] - s1[, 1] ^ 2))
  }

  # Cross-validation bandwidth selector (local linear)
  h.crossvalidation.ll <- function(n, x, y, hmin, hmax, hngrid) {
    crossvalue <- rep(0, hngrid)
    h <- seq(hmin, hmax, len = hngrid)
    for (j in 1:hngrid) {
      for (i in 1:n) {
        crossvalue[j] <- crossvalue[j] + (y[i] - locallinear(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^
          2
      }
    }
    crossvalue <- crossvalue / n
    #plot(h,crossvalue)
    h.crossvalidation.ll <- h[which.min(crossvalue)]
  }



  cat("Call:", "\n")
  print(match.call())
  DNAME <- deparse(substitute(c(x1, x2, y1, y2)))
  METHOD <- "A simple test for comparing regression curves versus one-sided alternatives"
  # DNAME <- "X"
  n1 <- length(x1)
  n2 <- length(x2)
  n <- n1 + n2
  kappa1 <- n1 / n
  kappa2 <- n2 / n

  # Functions p1 and p2
  p1.x1 <- sapply(x1, p)
  p2.x1 <- 1 - p1.x1
  p1.x2 <- sapply(x2, p)
  p2.x2 <- 1 - p1.x2

  # Functions w1 and w2
  w1.x1 <- sapply(x1, w1)
  w2.x1 <- sapply(x1, w2)
  w1.x2 <- sapply(x2, w1)
  w2.x2 <- sapply(x2, w2)

  h1.nw <- h.crossvalidation.nw(n1, x1, y1, 0.2, 0.4, 21)
  h2.nw <- h.crossvalidation.nw(n2, x2, y2, 0.2, 0.4, 21)

  m1.x1 <- nadarayawatson(n1, x1, y1, n1, x1, h1.nw)
  m1.x2 <- nadarayawatson(n1, x1, y1, n2, x2, h1.nw)
  m2.x1 <- nadarayawatson(n2, x2, y2, n1, x1, h2.nw)
  m2.x2 <- nadarayawatson(n2, x2, y2, n2, x2, h2.nw)

  m.x1 <- p1.x1 * m1.x1 + p2.x1 * m2.x1
  m.x2 <- p1.x2 * m1.x2 + p2.x2 * m2.x2

  sigma21.x1 <- nadarayawatson(n1, x1, y1 ^ 2, n1, x1, h1.nw) - m1.x1 ^ 2
  sigma21.x2 <- nadarayawatson(n1, x1, y1 ^ 2, n2, x2, h1.nw) - m1.x2 ^ 2
  sigma22.x1 <- nadarayawatson(n2, x2, y2 ^ 2, n1, x1, h2.nw) - m2.x1 ^ 2
  sigma22.x2 <- nadarayawatson(n2, x2, y2 ^ 2, n2, x2, h2.nw) - m2.x2 ^ 2

  # Densities
  h1.dens <- bw.ucv(x1, lower = 0.05, upper = 1.)
  h2.dens <- bw.ucv(x2, lower = 0.05, upper = 1.)

  f1.x1 <- kerneldensity(n1, x1, x1, h1.dens)
  f2.x1 <- kerneldensity(n2, x2, x1, h2.dens)
  f1.x2 <- kerneldensity(n1, x1, x2, h1.dens)
  f2.x2 <- kerneldensity(n2, x2, x2, h2.dens)

  # Test statistic

  eps01 <- y1 - m.x1
  eps02 <- y2 - m.x2

  f.x1 <- p1.x1 * w2.x1 * f2.x1 + p2.x1 * w1.x1 * f1.x1
  f.x2 <- p1.x2 * w2.x2 * f2.x2 + p2.x2 * w1.x2 * f1.x2

  sigma.asymp <- sqrt(kappa2 * mean(sigma21.x2 * f.x2 * p1.x2 * w2.x2 / f1.x2) +
                        kappa2 * mean(sigma21.x1 * f.x1 * p2.x1 * w1.x1 / f1.x1) +
                        kappa1 * mean(sigma22.x2 * f.x2 * p1.x2 * w2.x2 / f2.x2) +
                        kappa1 * mean(sigma22.x1 * f.x1 * p2.x1 * w1.x1 / f2.x1)
  )

  test.statistic <- sqrt(n1 * n2 / n) * (mean(eps02 * w2.x2) -
                                         mean(eps01 * w1.x1)) / sigma.asymp
  names( test.statistic) <- " test statistic"
  p.value <- 1 - pnorm(test.statistic)

  res <- list(p.value = p.value, data.name = DNAME,
              statistic =  test.statistic, method = METHOD)
  class(res) <- "htest"
  return(res)

}


#===============================================================================
# Example
#===============================================================================

n1 <- 100
n2 <- 200
x1 <- runif(n1)
y1 <- x1 + 0.25 * rnorm(n1)
x2 <- runif(n2)
y2 <- 0.2 + x2 + 0.5 * rnorm(n2)

comp2regr.onesided(x1, x2, y1, y2)


