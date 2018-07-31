#' @title function3
#' @description TODO
#' @param x1 TODO
#' @param y1 TODO
#' @param x2 TODO
#' @param y2 TODO
#' @param B bootstrap samples
#' @param bandwidths TODO
#' @param sigma.w TODO
#' @export
comp2condvar <- function(x1, y1, x2, y2, B = 1000, bandwidths = "cv", sigma.w = 1){

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
  # (xdata, ydata) and the bandwidth "h".
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
        crossvalue[j] <- crossvalue[j] + (y[i] - locallinear(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^ 2
      }
    }
    crossvalue <- crossvalue / n
    #plot(h, crossvalue)
    h.crossvalidation.ll <- h[which.min(crossvalue)]
  }

  # Functions related to the calculation of the test statistics based on the ecf

  # Weight function w: Normal(0,sigma)
  sigma.w <- 1

  w <- function(t, sigma.w) {
    dnorm(t, sd = sigma.w)
  }  # Normal

  # Function Iw (depends on the function w)
  Iw <- function(t, sigma.w) {
    exp(-0.5 * ((sigma.w * t) ^ 2))
  }  # w Normal

  # Function D2Iw (depends on the function w)
  D2Iw <- function(t, sigma.w) {
    ((sigma.w ^ 4) * t ^ 2 - sigma.w ^ 2) * exp(-0.5 * ((sigma.w * t) ^ 2))
  }  # w Normal

  # Estimator of the real and imaginary parts of the characteristic function
  C.hat <- function(t, n, x) {
    mean(cos(rep(t, n) * x))
  }

  S.hat <- function(t, n, x) {
    mean(sin(rep(t, n) * x))
  }

  # Test statistic
  Tn1.function <- function(n1, n2, n, eps1, eps01, eps2, eps02, sigma.w) {
      (sum(Iw(outer(eps1,  eps1,  "-"), sigma.w)) +
       sum(Iw(outer(eps01, eps01, "-"), sigma.w)) -
   2 * sum(Iw(outer(eps1,  eps01, "-"), sigma.w))) / n1 +
      (sum(Iw(outer(eps2,  eps2,  "-"), sigma.w)) +
       sum(Iw(outer(eps02, eps02, "-"), sigma.w)) -
   2 * sum(Iw(outer(eps2,  eps02, "-"), sigma.w))) / n2
  }


  Tn2.function <- function(n1, n2, n, eps1, eps01, eps2, eps02, sigma.w) {
    (sum(Iw(outer(c(eps1,  eps2),  c(eps1,  eps2),  "-"), sigma.w)) +
     sum(Iw(outer(c(eps01, eps02), c(eps01, eps02), "-"), sigma.w)) -
 2 * sum(Iw(outer(c(eps1,  eps2),  c(eps01, eps02), "-"), sigma.w))) / n
	}



  # Cramér-von Mises statistic
  emp.distr <- function(ndata, data, npoints, points) {
    sapply(data, function(data, points) {
      data <= points
    }, points = points) %*% rep(1, ndata) / ndata
  }

  cramer.vonmises <- function(n0, x0, n, x) {
    sum((emp.distr(n0, x0, n0, x0) - emp.distr(n, x, n0, x0)) ^ 2)
  }


#--------------------------------------------------------------------------------

  # Arguments
  # data: x1, y1, x2, y2
  # B (default = 1000)
  # h.m1, h.m2, h.sigma1, h.sigma2 (default = cross-validation)
  # sigma.w (default =1)

  bandwidths <- "cv"
  # Otherwise, specify h.m1, h.m2, h.sigma1, h.sigma2
  # Cross-validation specifications
  hmin   <- 0.05
  hmax   <- 0.5
  hngrid <- 46

  graphs <- 1
  bootstrap <- TRUE



  n1 <- length(x1)
  n2 <- length(x2)

  n  <- n1 + n2
  p1 <- n1 / n
  p2 <- n2 / n
  p.matrix <- diag(sqrt(c(p1, p2)))

  m0x1.hat <- rep(0, n1)
  m1x1.hat <- rep(0, n1)
  m2x1.hat <- rep(0, n1)
  m0x2.hat <- rep(0, n2)
  m1x2.hat <- rep(0, n2)
  m2x2.hat <- rep(0, n2)

  sigma0x1.hat <- rep(0, n1)
  sigma1x1.hat <- rep(0, n1)
  sigma1x2.hat <- rep(0, n2)
  sigma0x2.hat <- rep(0, n2)
  sigma2x2.hat <- rep(0, n2)
  sigma2x1.hat <- rep(0, n1)

  m0x1.hat.boot <- rep(0, n1)
  m1x1.hat.boot <- rep(0, n1)
  m2x1.hat.boot <- rep(0, n1)
  m0x2.hat.boot <- rep(0, n2)
  m1x2.hat.boot <- rep(0, n2)
  m2x2.hat.boot <- rep(0, n2)

  sigma0x1.hat.boot <- rep(0, n1)
  sigma1x1.hat.boot <- rep(0, n1)
  sigma1x2.hat.boot <- rep(0, n2)
  sigma0x2.hat.boot <- rep(0, n2)
  sigma2x2.hat.boot <- rep(0, n2)
  sigma2x1.hat.boot <- rep(0, n1)

  fmixx1.hat <- rep(0, n1)
  fmixx2.hat <- rep(0, n2)

  f1x1.hat <- rep(0, n1)
  f1x2.hat <- rep(0, n2)
  f2x1.hat <- rep(0, n1)
  f2x2.hat <- rep(0, n2)



  # Estimation of sigma1, sigma2 and sigma0

  if (bandwidths == "cv") {
    h.sigma1 <- h.crossvalidation.nw(n1, x1, y1, hmin, hmax, hngrid)
   #h.sigma1 <- npregbw(xdat = x1, ydat = y1 ^ 2, bwmethod = "cv.ls", lower = 0.05, upper = 0.5)$bw
    h.sigma2 <- h.crossvalidation.nw(n2, x2, y2, hmin, hmax, hngrid)
   #h.sigma2 <- npregbw(xdat = x2, ydat = y2 ^ 2, bwmethod = "cv.ls", lower = 0.05, upper = 0.5)$bw

    h.m1 <- h.crossvalidation.ll(n1, x1, y1, hmin, hmax, hngrid)
   #h.m1 <- npregbw(xdat = x1, ydat = y1, bwmethod = "cv.ls", lower = 0.05, upper = 0.5)$bw;
    h.m2 <- h.crossvalidation.ll(n2, x2, y2, hmin, hmax, hngrid)
   #h.m2 <- npregbw(xdat = x2, ydat = y2, bwmethod = "cv.ls", lower = 0.05, upper = 0.5)$bw;
  }

  sigma1x1.hat <- sqrt(abs(nadarayawatson(n1, x1, y1 ^ 2, n1, x1, h.sigma1) -
                           nadarayawatson(n1, x1, y1, n1, x1, h.sigma1) ^ 2))
  sigma1x2.hat <- sqrt(abs(nadarayawatson(n1, x1, y1 ^ 2, n2, x2, h.sigma1) -
                           nadarayawatson(n1, x1, y1, n2, x2,h.sigma1) ^ 2))
  sigma2x1.hat <- sqrt(abs(nadarayawatson(n2, x2, y2 ^ 2, n1, x1, h.sigma2) -
                           nadarayawatson(n2, x2, y2, n1, x1, h.sigma2) ^ 2))
  sigma2x2.hat <- sqrt(abs(nadarayawatson(n2, x2, y2 ^ 2, n2, x2, h.sigma2) -
                           nadarayawatson(n2, x2, y2, n2, x2, h.sigma2) ^ 2))

  pi1.x1 <- rep(p1, n1)
  pi2.x1 <- 1 - pi1.x1

  pi1.x2 <- rep(p1, n2)
  pi2.x2 <- 1 - pi1.x2

  sigma0x1.hat <- sqrt(pi1.x1 * (sigma1x1.hat ^ 2) + pi2.x1 * (sigma2x1.hat ^ 2))
  sigma0x2.hat <- sqrt(pi1.x2 * (sigma1x2.hat ^ 2) + pi2.x2 * (sigma2x2.hat ^ 2))

  # Estimation of m1 and m2

  m1x1.hat <- locallinear(n1, x1, y1, n1, x1, h.m1)
  m2x2.hat <- locallinear(n2, x2, y2, n2, x2, h.m2)

  # Estimation of the densities
  h.f1 <- stats::bw.ucv(x1, lower = 0.05, upper = 0.50)
  h.f2 <- stats::bw.ucv(x2, lower = 0.05, upper = 0.50)
  f1x1 <- kerneldensity(n1, x1, x1, h.f1)
  f1x2 <- kerneldensity(n1, x1, x2, h.f1)
  f2x1 <- kerneldensity(n2, x2, x1, h.f2)
  f2x2 <- kerneldensity(n2, x2, x2, h.f2)


  # Residuals
  eps1.hat  <- (y1 - m1x1.hat) / sigma1x1.hat
  eps01.hat <- (y1 - m1x1.hat) / sigma0x1.hat

  eps2.hat  <- (y2 - m2x2.hat) / sigma2x2.hat
  eps02.hat <- (y2 - m2x2.hat) / sigma0x2.hat


  # Test statistics
  Tn1  <- Tn1.function(n1, n2, n, eps1.hat, eps01.hat, eps2.hat, eps02.hat, sigma.w)
  Tn2  <- Tn2.function(n1, n2, n, eps1.hat, eps01.hat, eps2.hat, eps02.hat, sigma.w)

  KS.1 <- sqrt(n1) * stats::ks.test(eps01.hat, eps1.hat)$statistic + sqrt(n2) * stats::ks.test(eps02.hat, eps2.hat)$statistic
  KS.2 <- sqrt(n)  * stats::ks.test(c(eps01.hat, eps02.hat), c(eps1.hat, eps2.hat))$statistic

  CM.1 <- cramer.vonmises(n1, eps01.hat, n1, eps1.hat) + cramer.vonmises(n2, eps02.hat, n2, eps2.hat)
  CM.2 <- cramer.vonmises(n1 + n2, c(eps01.hat, eps02.hat), n1 + n2, c(eps1.hat, eps2.hat))


  # Estimation of the matrices involved in the approximation of the critical values of the test statistics Tn1 and CM.1

  # Standardization of residuals
  eps1.hat <- (eps1.hat - mean(eps1.hat)) / stats::sd(eps1.hat)
  eps2.hat <- (eps2.hat - mean(eps2.hat)) / stats::sd(eps2.hat)

  # Matrix A
  a1 <- -(sum((eps1.hat %*% t(eps1.hat)) * D2Iw(outer(eps1.hat,eps1.hat, "-"), sigma.w)) -
          eps1.hat%*%eps1.hat) / (4 * n1 * (n1 - 1))
  a2 <- -(sum((eps2.hat %*% t(eps2.hat)) * D2Iw(outer(eps2.hat,eps2.hat, "-"), sigma.w)) -
          eps2.hat%*%eps2.hat) / (4 * n2 * (n2 - 1))
  a.matrix <- diag(c(a1, a2))

  # Matrix F
  h.feps1 <- stats::bw.ucv(eps1.hat, lower = 0.10, upper = 1.00)
  h.feps2 <- stats::bw.ucv(eps2.hat, lower = 0.10, upper = 1.00)
  f1 <- 0.25 * mean((eps1.hat * kerneldensity(n1, eps1.hat, eps1.hat, h.feps1)) ^ 2)
  f2 <- 0.25 * mean((eps2.hat * kerneldensity(n2, eps2.hat, eps2.hat, h.feps2)) ^ 2)
  f.matrix <- diag(c(f1, f2))

  Eeps1 <- mean((eps1.hat ^ 2 - 1) ^ 2) / p1
  Eeps2 <- mean((eps2.hat ^ 2 - 1) ^ 2) / p2

  sigma11 <- p1 * (Eeps1 * mean((pi1.x1 - 1) ^ 2) + Eeps2 * mean((pi2.x2 * f1x2 / f2x2) ^ 2))
  sigma22 <- p2 * (Eeps1 * mean((pi1.x1 * f2x1 / f1x1) ^ 2) + Eeps2 * mean((pi2.x2 - 1) ^ 2) )
  sigma12 <- sqrt(p1 * p2)*(Eeps1 * mean((pi1.x1 - 1) * pi1.x1 * f2x1 / f1x1) + Eeps2 * mean((pi2.x2 * f1x2 / f2x2) * (pi2.x2 - 1)))
  sigma.matrix <- matrix(c(sigma11, sigma12, sigma12, sigma22), nrow =  2)

  # Eigenvalues: coefficients beta and tseta
  beta  <- eigen(a.matrix %*% sigma.matrix)$values
  beta1 <- beta[1]
  beta2 <- beta[2]

  tseta  <- eigen(f.matrix %*% sigma.matrix)$values
  tseta1 <- tseta[1]
  tseta2 <- tseta[2]

  # p-values of Tn1 and CM.1 based on the asymptotic null distribution
  pvalue.Tn1.asym  <- mean(Tn1 < beta[1]   * stats::rchisq(10 ^ 6, 1) + beta[2]  * stats::rchisq(10 ^ 6, 1))
  pvalue.CM.1.asym <- mean(CM.1 < tseta[1] * stats::rchisq(10 ^ 6, 1) + tseta[2] * stats::rchisq(10 ^ 6, 1))

  # BOOTSTAP
  if (bootstrap == TRUE){

    KS.1.boot <- rep(0., B)
    KS.2.boot <- rep(0., B)

    CM.1.boot <- rep(0., B)
    CM.2.boot <- rep(0., B)

    Tn1.boot <- rep(0, B)
    Tn2.boot <- rep(0, B)

    a1smooth <- 0.
    a2smooth <- 0.

    # standardized residuals
    eps1.stand <- (eps1.hat - mean(eps1.hat)) / stats::sd(eps1.hat)
    eps2.stand <- (eps2.hat - mean(eps2.hat)) / stats::sd(eps2.hat)

    eps1.boot.matrix <- matrix(sqrt(1 - a1smooth ^ 2) * sample(eps1.hat, size = B * n1, replace = TRUE) +
                               a1smooth * stats::rnorm(B * n1), nrow = B, ncol = n1)
    eps2.boot.matrix <- matrix(sqrt(1 - a2smooth ^ 2) * sample(eps2.hat, size = B * n2, replace = TRUE) +
                               a2smooth * stats::rnorm(B * n2), nrow = B, ncol = n2)

    for (ib in 1:B) {

  #	print(c("ib",ib))

      y1.boot <- m0x1.hat + sigma0x1.hat * eps1.boot.matrix[ib, ]
      y2.boot <- m0x2.hat + sigma0x2.hat * eps2.boot.matrix[ib, ]

      sigma1x1.hat.boot <- sqrt(abs(nadarayawatson(n1, x1, y1.boot ^ 2, n1, x1, h.sigma1) -
                                    nadarayawatson(n1, x1, y1.boot, n1, x1, h.sigma1) ^ 2))
      sigma1x2.hat.boot <- sqrt(abs(nadarayawatson(n1, x1, y1.boot ^ 2, n2, x2, h.sigma1) -
                                    nadarayawatson(n1, x1, y1.boot, n2, x2, h.sigma1) ^ 2))

      sigma2x1.hat.boot <- sqrt(abs(nadarayawatson(n2,x2, y2.boot ^ 2, n1, x1, h.sigma2) -
                                    nadarayawatson(n2, x2, y2.boot, n1, x1, h.sigma2) ^ 2))
      sigma2x2.hat.boot <- sqrt(abs(nadarayawatson(n2, x2, y2.boot ^ 2, n2, x2, h.sigma2) -
                                    nadarayawatson(n2, x2, y2.boot, n2, x2, h.sigma2) ^ 2))

      sigma0x1.hat.boot <- sqrt(pi1.x1 * (sigma1x1.hat.boot ^ 2) + pi2.x1 * (sigma2x1.hat.boot ^ 2))
      sigma0x2.hat.boot <- sqrt(pi1.x2 * (sigma1x2.hat.boot ^ 2) + pi2.x2 * (sigma2x2.hat.boot ^ 2))

      m1x1.hat.boot <- locallinear(n1, x1, y1.boot, n1, x1, h.m1)
      m2x2.hat.boot <- locallinear(n2, x2, y2.boot, n2, x2, h.m2)

      # Bootstrap residuals
      eps1.hat.boot  <- (y1.boot - m1x1.hat.boot) / sigma1x1.hat.boot
      eps01.hat.boot <- (y1.boot - m1x1.hat.boot) / sigma0x1.hat.boot

      eps2.hat.boot  <- (y2.boot - m2x2.hat.boot) / sigma2x2.hat.boot
      eps02.hat.boot <- (y2.boot - m2x2.hat.boot) / sigma0x2.hat.boot

      # Test statistics
      Tn1.boot[ib]  <- Tn1.function(n1, n2, n, eps1.hat.boot, eps01.hat.boot,
                                    eps2.hat.boot, eps02.hat.boot, sigma.w)
      Tn2.boot[ib]  <- Tn2.function(n1, n2, n, eps1.hat.boot, eps01.hat.boot
                                    , eps2.hat.boot, eps02.hat.boot, sigma.w)

      KS.1.boot[ib] <- sqrt(n1) * stats::ks.test(eps01.hat.boot, eps1.hat.boot)$statistic +
                       sqrt(n2) * stats::ks.test(eps02.hat.boot, eps2.hat.boot)$statistic
      KS.2.boot[ib] <- sqrt(n)  * stats::ks.test(c(eps01.hat.boot, eps02.hat.boot),
                                                c(eps1.hat.boot, eps2.hat.boot))$statistic

      CM.1.boot[ib] <- cramer.vonmises(n1, eps01.hat.boot, n1, eps1.hat.boot) +
                       cramer.vonmises(n2, eps02.hat.boot, n2, eps2.hat.boot)
      CM.2.boot[ib] <- cramer.vonmises(n1 + n2, c(eps01.hat.boot, eps02.hat.boot),
                                       n1 + n2, c(eps1.hat.boot, eps2.hat.boot))

  }  #ib

  # Bootstrap p-values
  pvalue.Tn1.boot  <- mean(Tn1 < Tn1.boot)
  pvalue.Tn2.boot  <- mean(Tn2 < Tn2.boot)

  pvalue.KS.1.boot <- mean(KS.1 < KS.1.boot)
  pvalue.KS.2.boot <- mean(KS.2 < KS.2.boot)

  pvalue.CM.1.boot <- mean(CM.1 < CM.1.boot)
  pvalue.CM.2.boot <- mean(CM.2 < CM.2.boot)

} # bootstrap ?


  # Output:
  print(list(pvalue.Tn1.boot  = pvalue.Tn1.boot,  pvalue.Tn1.asym  = pvalue.Tn1.asym,
            pvalue.Tn2.boot  = pvalue.Tn2.boot,  pvalue.KS.1.boot = pvalue.KS.1.boot,
            pvalue.KS.2.boot = pvalue.KS.2.boot, pvalue.CM.1.boot = pvalue.CM.1.boot,
            pvalue.CM.1.asym = pvalue.CM.1.asym, pvalue.CM.2.boot = pvalue.CM.2.boot))


  r <- list(pvalue.Tn1.boot  = pvalue.Tn1.boot,  pvalue.Tn1.asym  = pvalue.Tn1.asym,
  pvalue.Tn2.boot  = pvalue.Tn2.boot,  pvalue.KS.1.boot = pvalue.KS.1.boot,
  pvalue.KS.2.boot = pvalue.KS.2.boot, pvalue.CM.1.boot = pvalue.CM.1.boot,
  pvalue.CM.1.asym = pvalue.CM.1.asym, pvalue.CM.2.boot = pvalue.CM.2.boot, x1 = x1,
  x2 = x2, y1 = y1, y2 = y2, m1x1.hat = m1x1.hat, m0x1.hat = m0x1.hat,
  m2x2.hat = m2x2.hat, m0x2.hat = m0x2.hat, sigma1x1.hat = sigma1x1.hat,
  sigma2x2.hat = sigma2x2.hat, eps1.hat = eps1.hat, eps01.hat = eps01.hat,
  eps2.hat = eps2.hat, eps02.hat = eps02.hat)

  class(r) <- c('list', 'comp2condvar')
  return(r)

  # pvalue.Tn1.boot
  # pvalue.Tn1.asym
  # pvalue.Tn2.boot
  # pvalue.KS.1.boot
  # pvalue.KS.2.boot
  # pvalue.CM.1.boot
  # pvalue.CM.1.asym
  # pvalue.CM.2.boot
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
# res <- comp2condvar(x1, y1, x2, y2, B = 100)


#===============================================================================
# Plot function
#===============================================================================

# if (graphs == 1){
#
#   par(mfrow = c(2, 2))
#
#   plot(c(x1, x2),c(y1, y2), type = "n", xlab = "covariate", ylab = "response",
#        main = "resgression functions")
#   points(x1, y1, col = "red")
#   points(x2, y2, col = "blue")
#   lines(sort(x1), m1x1.hat[order(x1)], col = "red")
#   lines(sort(x2), m2x2.hat[order(x2)], col = "blue")
#
#   plot(c(x1, x2), c(y1, y2), type = "n", xlab = "covariate",
#      ylab = "sigma", main = "conditional variance functions")
#   lines(sort(x1), sigma1x1.hat[order(x1)], col = "red")
#   lines(sort(x1), sigma0x1.hat[order(x1)], col = "red", lty = 2)
#   lines(sort(x2), sigma2x2.hat[order(x2)], col = "blue")
#   lines(sort(x2), sigma0x2.hat[order(x2)], col = "blue", lty = 2)
#
#   plot(ecdf(eps1.hat), col = "red", main = "F_eps1")
#   plot(ecdf(eps01.hat), col = "black", add = T)
#
#   plot(ecdf(eps2.hat), col = "blue", main = "F_eps2")
#   plot(ecdf(eps02.hat), col = "black", add = T)
#
# }   #graphs




