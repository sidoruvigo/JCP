#' @title function2
#' @description description of the function
#' @param x1 TODO
#' @param y1 TODO
#' @param x2 TODO
#' @param y2 TODO
#' @param B bootstrap samples
#' @param bandwidths TODO
#' @param sigma.w TODO
#' @importFrom graphics lines par plot points
#' @importFrom stats dnorm ecdf pnorm
#' @export
comp2regr.ecdf <- function(x1, y1, x2, y2, B = 1000, bandwidths = "cv", sigma.w = 1) {

  # GENERAL FUNCTIONS (kernel-based estimators)

  # Kernel Epanechnikov (density)
  kernel <- function(u) {
    (0.75 * (1 - u ^ 2)) * (u < 1) * (u > -1)
  }

  # Kernel density estimator. This function calculates the kernel density
  # estimator on the collection of "npoints" points "points" based on the "ndata"
  # observations xdata and the bandwidth "h".
  kerneldensity <- function(ndata, data, points, h) {
    rowSums(kernel(outer(points, data, "-") / h)) / (ndata * h)
  }

  # Nadaraya-Watson estimator. This function calculates the N-W estimator on the
  # collection of "npoints" points "points" based on the "ndata" observations
  # (xdata,ydata) and the bandwidth "h".
  nadarayawatson <- function(ndata, xdata, ydata, npoints, points, h) {
    as.vector({
      matk <- kernel((points %*% t(rep(1, ndata)) - t(xdata %*% t(rep(1, npoints)))) / h)
     (matk %*% ydata) / (matk %*% rep(1, ndata))
    })
  }

  # Cross-validation bandwidth selector (Nadaraya-Watson)
  h.crossvalidation.nw <- function(n, x, y, hmin, hmax, hngrid) {
    crossvalue <- rep(0, hngrid)
    h <- seq(hmin, hmax, len = hngrid)
    for (j in 1:hngrid) {
      for (i in 1:n) {
        crossvalue[j] <- crossvalue[j] + (y[i] - nadarayawatson(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^ 2
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
        crossvalue[j] <- crossvalue[j] + (y[i] - locallinear(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^ 2
      }
    }
    crossvalue <- crossvalue / n
    #plot(h,crossvalue)
    h.crossvalidation.ll <- h[which.min(crossvalue)]
  }


  # Cramér-von Mises statistic
  emp.distr <- function(ndata, data, npoints, points) {
    sapply(data, function(data, points) {
      data <= points
    }, points = points) %*% rep(1, ndata) / ndata
  }

  cramer.vonmises <- function(nx, x, ny, y) {
    sum((emp.distr(nx, x, nx + ny, c(x, y)) - emp.distr(ny, y, nx + ny, c(x, y))) ^ 2)
  }


#-------------------------------------------------------------------------------

  # Arguments
  # data: x1, y1, x2, y2
  # B (default=1000)
  # h.m1, h.m2, h.sigma1, h.sigma2 (default = cross-validation)

  # bandwidths <- "cv"

  # Otherwise, specify h.m1, h.m2, h.sigma1, h.sigma2

  # Cross-validation specifications
  hmin   <- 0.05
  hmax   <- 0.5
  hngrid <- 46


  # graphs <- 1
  # B <- 200


  ########################
  ########################


  n1 <- length(x1)
  n2 <- length(x2)
  n  <- n1 + n2

  p1 <- n1 / n
  p2 <- n2 / n

  m0x1.hat <- rep(0, n1)
  m1x1.hat <- rep(0, n1)
  m2x1.hat <- rep(0, n1)
  m0x2.hat <- rep(0, n2)
  m1x2.hat <- rep(0, n2)
  m2x2.hat <- rep(0, n2)

  sigma1x1.hat <- rep(0, n1)
  sigma1x2.hat <- rep(0, n2)
  sigma2x2.hat <- rep(0, n2)
  sigma2x1.hat <- rep(0, n1)

  m0x1.hat.boot <- rep(0, n1)
  m1x1.hat.boot <- rep(0, n1)
  m2x1.hat.boot <- rep(0, n1)
  m0x2.hat.boot <- rep(0, n2)
  m1x2.hat.boot <- rep(0, n2)
  m2x2.hat.boot <- rep(0, n2)

  sigma1x1.hat.boot <- rep(0, n1)
  sigma1x2.hat.boot <- rep(0, n2)
  sigma2x2.hat.boot <- rep(0, n2)
  sigma2x1.hat.boot <- rep(0, n1)

  fmixx1.hat <- rep(0, n1)
  fmixx2.hat <- rep(0, n2)

  f1x1.hat <- rep(0, n1)
  f1x2.hat <- rep(0, n2)
  f2x1.hat <- rep(0, n1)
  f2x2.hat <- rep(0, n2)

  if (bandwidths == "cv") {

    h.sigma1 <- h.crossvalidation.nw(n1, x1, y1, hmin, hmax, hngrid)
   #h.sigma1 <- npregbw(xdat = x1, ydat = y1, bwmethod = "cv.ls", kernel = "epanech", regtype = "lc")$bw
    h.sigma2 <- h.crossvalidation.nw(n2, x2, y2, hmin, hmax, hngrid)
   #h.sigma2 <- npregbw(xdat = x2, ydat = y2, bwmethod = "cv.ls", kernel = "epanech", regtype = "lc")$bw

    h.m1 <- h.crossvalidation.ll(n1, x1, y1, hmin, hmax, hngrid)
   #h.m1 <- npregbw(xdat = x1, ydat = y1, bwmethod = "cv.ls", kernel = "epanech", regtype = "ll")$bw
    h.m2 <- h.crossvalidation.ll(n2, x2, y2, hmin, hmax, hngrid)
   #h.m1 <- npregbw(xdat = x2, ydat = y2, bwmethod = "cv.ls", kernel = "epanech", regtype = "ll")$bw

  } else {

    h.sigma1 <- bandwidths[1]
    h.sigma2 <- bandwidths[2]

    h.m1 <- bandwidths[3]
    h.m2 <- bandwidths[4]

  }

  # Estimation of m1, m2, m0, sigma1, sigma2

  sigma1x1.hat <- sqrt(abs(nadarayawatson(n1, x1, y1 ^ 2, n1, x1, h.sigma1) -
                           nadarayawatson(n1, x1, y1, n1, x1, h.sigma1) ^ 2))
  sigma2x2.hat <- sqrt(abs(nadarayawatson(n2, x2, y2 ^ 2, n2, x2, h.sigma2) -
                           nadarayawatson(n2, x2, y2, n2, x2, h.sigma2) ^ 2))

  m1x1.hat <- locallinear(n1, x1, y1, n1, x1, h.m1)
  m1x2.hat <- locallinear(n1, x1, y1, n2, x2, h.m1)
  m2x2.hat <- locallinear(n2, x2, y2, n2, x2, h.m2)
  m2x1.hat <- locallinear(n2, x2, y2, n1, x1, h.m2)

  h.f1 <- stats::bw.nrd0(x1)
  f1x1 <- kerneldensity(n1, x1, x1, h.f1)
  f1x2 <- kerneldensity(n1, x1, x2, h.f1)

  h.f2 <- stats::bw.nrd0(x2)
  f2x1 <- kerneldensity(n2, x2, x1, h.f2)
  f2x2 <- kerneldensity(n2, x2, x2, h.f2)

  h.fmix <- stats::bw.nrd0(c(x1, x2))
  fmixx1 <- kerneldensity(n1 + n2, c(x1, x2), x1, h.fmix)
  fmixx2 <- kerneldensity(n1 + n2, c(x1, x2), x2, h.fmix)

  m0x1.hat <- (p1 * f1x1 * m1x1.hat + p2 * f2x1 * m2x1.hat) / fmixx1
  m0x2.hat <- (p1 * f1x2 * m1x2.hat + p2 * f2x2 * m2x2.hat) / fmixx2


  # Residuals
  eps1.hat  <- (y1 - m1x1.hat) / sigma1x1.hat
  eps01.hat <- (y1 - m0x1.hat) / sigma1x1.hat

  eps2.hat  <- (y2 - m2x2.hat) / sigma2x2.hat
  eps02.hat <- (y2 - m0x2.hat) / sigma2x2.hat

  # Test statistics
  KS.1 <- sqrt(n1) * stats::ks.test(eps01.hat, eps1.hat)$statistic + sqrt(n2) * stats::ks.test(eps02.hat, eps2.hat)$statistic
  KS.2 <- sqrt(n)  * stats::ks.test(c(eps01.hat, eps02.hat), c(eps1.hat, eps2.hat))$statistic

  CM.1 <- cramer.vonmises(n1, eps01.hat, n1, eps1.hat) + cramer.vonmises(n2, eps02.hat, n2, eps2.hat)
  CM.2 <- cramer.vonmises(n1 + n2, c(eps01.hat, eps02.hat), n1 + n2, c(eps1.hat, eps2.hat))

  # p-values based on BOOTSTAP (smooth bootstrap of residuals. Statistica Sinica 2007)
  KS.1.boot <- rep(0., B)
  KS.2.boot <- rep(0., B)
  CM.1.boot <- rep(0., B)
  CM.2.boot <- rep(0., B)

  a1smooth <- 0.
  a2smooth <- 0.

  # Standardized residuals
  eps1.stand <- (eps1.hat - mean(eps1.hat)) / stats::sd(eps1.hat)
  eps2.stand <- (eps2.hat - mean(eps2.hat)) / stats::sd(eps2.hat)

  eps1.boot.matrix <- matrix(sqrt(1 - a1smooth ^ 2) *
                             sample(eps1.stand, size = B * n1, replace = TRUE) +
                             a1smooth * stats::rnorm(B * n1), nrow = B, ncol = n1)
  eps2.boot.matrix <- matrix(sqrt(1 - a2smooth ^ 2) *
                             sample(eps2.stand, size = B * n2, replace = TRUE) +
                             a2smooth * stats::rnorm(B * n2), nrow = B, ncol = n2)

  for (ib in 1:B){

  #	print(c("ib",ib))

    y1.boot <- m0x1.hat + sigma1x1.hat * eps1.boot.matrix[ib, ]
    y2.boot <- m0x2.hat + sigma2x2.hat * eps2.boot.matrix[ib, ]

    sigma1x1.hat.boot <- sqrt(abs(nadarayawatson(n1, x1, y1.boot ^ 2, n1, x1, h.sigma1) -
                                  nadarayawatson(n1, x1, y1.boot, n1, x1, h.sigma1) ^ 2))
    sigma2x2.hat.boot <- sqrt(abs(nadarayawatson(n2, x2, y2.boot ^ 2, n2, x2, h.sigma2) -
                                  nadarayawatson(n2, x2, y2.boot, n2, x2, h.sigma2) ^ 2))


    m1x1.hat.boot <- locallinear(n1, x1, y1.boot, n1, x1, h.m1)
    m1x2.hat.boot <- locallinear(n1, x1, y1.boot, n2, x2, h.m1)
    m2x2.hat.boot <- locallinear(n2, x2, y2.boot, n2, x2, h.m2)
    m2x1.hat.boot <- locallinear(n2, x2, y2.boot, n1, x1, h.m2)

    m0x1.hat.boot <- (p1 * f1x1 * m1x1.hat.boot + p2 * f2x1 * m2x1.hat.boot) / fmixx1
    m0x2.hat.boot <- (p1 * f1x2 * m1x2.hat.boot + p2 * f2x2 * m2x2.hat.boot) / fmixx2


  # Bootstrap residuals
    eps1.hat.boot  <- (y1.boot - m1x1.hat.boot) / sigma1x1.hat.boot
    eps01.hat.boot <- (y1.boot - m0x1.hat.boot) / sigma1x1.hat.boot

    eps2.hat.boot  <- (y2.boot - m2x2.hat.boot) / sigma2x2.hat.boot
    eps02.hat.boot <- (y2.boot - m0x2.hat.boot) / sigma2x2.hat.boot

  # Test statistics
    KS.1.boot[ib] <- sqrt(n1) * stats::ks.test(eps01.hat.boot, eps1.hat.boot)$statistic +
                     sqrt(n2) * stats::ks.test(eps02.hat.boot, eps2.hat.boot)$statistic
    KS.2.boot[ib] <- sqrt(n) * stats::ks.test(c(eps01.hat.boot, eps02.hat.boot),
                                              c(eps1.hat.boot, eps2.hat.boot))$statistic

    CM.1.boot[ib] <- cramer.vonmises(n1, eps01.hat.boot, n1, eps1.hat.boot) +
                     cramer.vonmises(n2, eps02.hat.boot, n2, eps2.hat.boot)
    CM.2.boot[ib] <- cramer.vonmises(n1 + n2, c(eps01.hat.boot, eps02.hat.boot), n1 +
                                     n2, c(eps1.hat.boot,eps2.hat.boot))

  }  #ib

  pvalue.KS.1.boot <- mean(KS.1 < KS.1.boot)
  pvalue.CM.1.boot <- mean(CM.1 < CM.1.boot)
  pvalue.KS.2.boot <- mean(KS.2 < KS.2.boot)
  pvalue.CM.2.boot <- mean(CM.2 < CM.2.boot)

  # Output:

  # pvalue.KS.1.boot
  # pvalue.CM.1.boot
  # pvalue.KS.2.boot
  # pvalue.CM.2.boot

  print(list(pvalue.KS.1.boot = pvalue.KS.1.boot, pvalue.CM.1.boot = pvalue.CM.1.boot,
             pvalue.KS.2.boot = pvalue.KS.2.boot, pvalue.CM.2.boot = pvalue.CM.2.boot))

  r <- list(pvalue.KS.1.boot = pvalue.KS.1.boot, pvalue.CM.1.boot = pvalue.CM.1.boot,
             pvalue.KS.2.boot = pvalue.KS.2.boot, pvalue.CM.2.boot = pvalue.CM.2.boot, x1 = x1,
             x2 = x2, y1 = y1, y2 = y2, m1x1.hat = m1x1.hat, m0x1.hat = m0x1.hat,
             m2x2.hat = m2x2.hat, m0x2.hat = m0x2.hat, sigma1x1.hat = sigma1x1.hat,
             sigma2x2.hat = sigma2x2.hat, eps1.hat = eps1.hat, eps01.hat = eps01.hat,
             eps2.hat = eps2.hat, eps02.hat = eps02.hat)

  class(r) <- c('list', 'comp2regr.ecdf')
  return(r)
}




#===============================================================================
# Example
#===============================================================================

# n1 <- 100
# n2 <- 200
# x1 <- runif(n1)
# y1 <- x1 + 0.25 * rnorm(n1)
# x2 <- runif(n2)
# y2 <- x2 + 0.25 * rnorm(n2)
#
# res <- comp2regr.ecdf(x1, y1, x2, y2)


# Plots (optional)
# if (graphs == 1) {
  # par(mfrow = c(2, 2))
  #
  # plot(c(x1, x2), c(y1, y2), type = "n", xlab = "covariate",
  #   ylab = "response", main = "resgression functions")
  # points(x1, y1, col = "red")
  # points(x2, y2, col = "blue")
  # lines(sort(x1), m1x1.hat[order(x1)], col = "red")
  # lines(sort(x1), m0x1.hat[order(x1)], col = "red", lty = 2)
  # lines(sort(x2), m2x2.hat[order(x2)], col = "blue")
  # lines(sort(x2), m0x2.hat[order(x2)], col = "blue", lty = 2)
  #
  # plot(c(x1, x2), c(y1, y2), type = "n", xlab = "covariate",
  #   ylab = "sigma", main = "conditional variance functions")
  # lines(sort(x1), sigma1x1.hat[order(x1)], col = "red")
  # lines(sort(x2), sigma2x2.hat[order(x2)], col = "blue")
  #
  # plot(ecdf(eps1.hat), col = "red", main = "F_eps1")
  # plot(ecdf(eps01.hat), col = "black", add = T)
  #
  # plot(ecdf(eps2.hat), col = "blue", main = "F_eps2")
  # plot(ecdf(eps02.hat), col = "black", add = T)

# }   #graphs

