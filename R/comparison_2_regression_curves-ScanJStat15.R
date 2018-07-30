#' @title function1
#' @export
function1 <- function(x1, y1, x2, y2, B = 100, bandwidths = "cv", sigma.w = 1){

  if(is.null(B)) {
    B <- 1000
  }

  if (is.null(bandwidths)) {
    bandwidths <- "cv"
  }


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
    as.vector({matk = kernel((points %*% t(rep(1, ndata)) - t(xdata %*% t(rep(1, npoints)))) / h)
    (matk %*% ydata) / (matk %*% rep(1, ndata))
    })
  }

  # Cross-validation bandwidth selector (Nadaraya-Watson)
  h.crossvalidation.nw <- function(n, x, y, hmin, hmax, hngrid) {
    crossvalue = rep(0, hngrid)
    h = seq(hmin, hmax, len = hngrid)
    for (j in 1:hngrid) {
      for (i in 1:n) {
        crossvalue[j] = crossvalue[j] + (y[i] - nadarayawatson(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^
          2
      }
    }
    crossvalue = crossvalue / n
    #plot(h,crossvalue)
    h.crossvalidation.nw = h[which.min(crossvalue)]
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




  # Functions related to the calculation of the test statistic Tn1 and Tn2

  # weight function w: Normal(0,sigma)
  # sigma.w <- 1
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

  # Test statistics T1 and T2
  Tn1.function <- function(n1, n2, n, eps1, eps01, eps2, eps02, sigma.w) {
    (sum(Iw(outer(eps1, eps1, "-"), sigma.w)) + sum(Iw(outer(eps01, eps01, "-"), sigma.w)) -
       2 * sum(Iw(outer(eps1, eps01, "-"), sigma.w))) / n1 + (sum(Iw(outer(eps2, eps2, "-"), sigma.w)) +
                                                                sum(Iw(outer(eps02, eps02, "-"), sigma.w)) - 2 * sum(Iw(outer(eps2, eps02, "-"), sigma.w))) /n2
  }

  Tn2.function <- function(n1, n2, n, eps1, eps01, eps2, eps02, sigma.w) {
    (sum(Iw(outer(c(eps1, eps2), c(eps1, eps2), "-"), sigma.w)) +
       sum(Iw(outer(c(eps01, eps02), c(eps01, eps02), "-"), sigma.w)) -
       2 * sum(Iw(outer(c(eps1, eps2), c(eps01, eps02), "-"), sigma.w))) / n
  }


  # Arguments
  # data: x1, y1, x2, y2
  # B (default=1000)
  # h.m1, h.m2, h.sigma1, h.sigma2 (default=cross-validation)
  # sigma.w (default =1)

  # bandwidths = "cv"
  # Otherwise, specify h.m1, h.m2, h.sigma1, h.sigma2
  # Cross-validation specifications
  hmin   <- 0.05
  hmax   <- 0.5
  hngrid <- 46

  # graphs <- 1
  bootstrap <- TRUE
  # B <- 200

  ########################
  ########################

  n1 <- length(x1)
  n2 <- length(x2)

  n <- n1 + n2

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
    #h.sigma1 = npregbw(xdat=x1,ydat=y1,bwmethod="cv.ls",kernel="epanech",regtype="lc")$bw
    h.sigma2 <- h.crossvalidation.nw(n2, x2, y2, hmin, hmax, hngrid)
    #h.m2.nw = npregbw(xdat=x2,ydat=y2,bwmethod="cv.ls",kernel="epanech",regtype="lc")$bw

    h.m1 <- h.crossvalidation.ll(n1, x1, y1, hmin, hmax, hngrid)
    #h.m1 = npregbw(xdat=x1,ydat=y1,bwmethod="cv.ls",kernel="epanech",regtype="ll")$bw
    h.m2 <- h.crossvalidation.ll(n2, x2, y2, hmin, hmax, hngrid)
    #h.m1 = npregbw(xdat=x2,ydat=y2,bwmethod="cv.ls",kernel="epanech",regtype="ll")$bw

  }

  # Estimation of m1, m2, m0, sigma1, sigma2

  sigma1x1.hat <- sqrt(abs(nadarayawatson(n1, x1, y1 ^ 2, n1, x1, h.sigma1) -
                           nadarayawatson(n1, x1, y1, n1, x1, h.sigma1) ^ 2))
  sigma1x2.hat <- sqrt(abs(nadarayawatson(n1, x1, y1 ^ 2, n2, x2, h.sigma1) -
                           nadarayawatson(n1, x1, y1, n2, x2, h.sigma1)^2))

  sigma2x1.hat <- sqrt(abs(nadarayawatson(n2, x2, y2 ^ 2, n1, x1, h.sigma2) -
                           nadarayawatson(n2, x2, y2, n1, x1, h.sigma2) ^ 2))
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

  Tn1 <- Tn1.function(n1, n2, n, eps1.hat, eps01.hat, eps2.hat, eps02.hat, sigma.w)

  Tn2 <- Tn2.function(n1, n2, n, eps1.hat, eps01.hat, eps2.hat, eps02.hat, sigma.w)

  # p-value based on the asymptotic null distribution of Tn1 (Scand. J. Stat, 2015)

  # Estimation of the matrix A
  a.matrix <- diag(c(-(sum(D2Iw(outer(eps1.hat, eps1.hat, "-"), sigma.w)) - n1) / (n1 * (n1 - 1)), -(sum(D2Iw(outer(eps2.hat, eps2.hat, "-"), sigma.w)) - n2) / (n2 * (n2 - 1))))

  # Estimation of matrix Sigma
  sigma11 <- 1 - 2 * p1 * mean(f1x1 / fmixx1) + p1 * (p1 * mean((f1x1/fmixx1)^2) + p2*mean((f1x2*sigma2x2.hat/(fmixx2*sigma1x2.hat))^2) )
  sigma22 <- 1 - 2 * p2 * mean(f2x2 / fmixx2) + p2 * (p1 * mean((f2x1*sigma1x1.hat/(fmixx1*sigma2x1.hat))^2) + p2*mean((f2x2/fmixx2)^2) )

  sigma12 <- sqrt(p1 * p2) * (p1 * mean(f1x1 * f2x1 * sigma1x1.hat / (sigma2x1.hat * (fmixx1 ^ 2))) +
                   p2 * mean(f1x2 * f2x2 * sigma2x2.hat / (sigma1x2.hat * (fmixx2 ^ 2))) -
                     mean(sigma2x2.hat * f1x2 / (sigma1x2.hat * fmixx2)) -
                     mean(sigma1x1.hat * f2x1 / (sigma2x1.hat * fmixx1)))

  sigma.matrix <- matrix(c(sigma11, sigma12, sigma12, sigma22), nrow = 2)

  # Estimation of betas
  beta  <- eigen(a.matrix %*% sigma.matrix)$values
  beta1 <- beta[1]
  beta2 <- beta[2]

  # p-value
  pvalue.Tn1.asym <- mean(Tn1 < beta[1] * rchisq(10 ^ 6, 1) + beta[2] * stats::rchisq(10 ^ 6, 1))



  # p-values based on BOOTSTAP (smooth bootstrap of residuals. Statistica Sinica 2007, Scand. J. Stat. 2015)
  if (bootstrap == TRUE){

    Tn2.boot <- rep(0., B)

    a1smooth <- 0.
    a2smooth <- 0.

    # Standardized residuals

    eps1.stand <- (eps1.hat - mean(eps1.hat)) / stats::sd(eps1.hat)
    eps2.stand <- (eps2.hat - mean(eps2.hat)) / stats::sd(eps2.hat)

    eps1.boot.matrix <- matrix(sqrt(1 - a1smooth ^ 2) *
                               sample(eps1.stand, size = B * n1, replace = TRUE) +
                               a1smooth * stats::rnorm(B * n1), nrow = B,  ncol = n1)
    eps2.boot.matrix <- matrix(sqrt(1 - a2smooth ^ 2) *
                               sample(eps2.stand, size = B * n2, replace = TRUE) +
                               a2smooth * stats::rnorm(B * n2), nrow = B, ncol = n2)

  for (ib in 1:B){

  #	print(c("ib",ib))

    y1.boot <- m0x1.hat + sigma1x1.hat * eps1.boot.matrix[ib, ]
	  y2.boot <- m0x2.hat + sigma2x2.hat * eps2.boot.matrix[ib, ]

	  sigma1x1.hat.boot <- sqrt(abs(nadarayawatson(n1, x1, y1.boot ^ 2, n1, x1, h.sigma1) -
	                                nadarayawatson(n1, x1, y1.boot, n1, x1, h.sigma1) ^ 2))
	  sigma2x2.hat.boot <- sqrt(abs(nadarayawatson(n2, x2, y2.boot ^ 2, n2, x2, h.sigma1) -
	                                nadarayawatson(n2, x2, y2.boot, n2, x2, h.sigma2)^2))

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

  Tn1.boot[ib] <- Tn1.function(n1, n2, n, eps1.hat.boot, eps01.hat.boot,
                               eps2.hat.boot, eps02.hat.boot, sigma.w)
  Tn2.boot[ib] <- Tn2.function(n1, n2, n, eps1.hat.boot, eps01.hat.boot,
                               eps2.hat.boot, eps02.hat.boot, sigma.w)

  }  #ib

    pvalue.Tn1.boot <- mean(Tn1 < Tn1.boot)
    pvalue.Tn2.boot <- mean(Tn2 < Tn2.boot)

  }   # boostrap?


  # Output:

  r1 <- list(Tn1.asym = pvalue.Tn1.asym,
             Tn1.boot = pvalue.Tn1.boot,
             Tn2.boot = pvalue.Tn2.boot)
  print(r1)

  r <- list(Tn1.asym = pvalue.Tn1.asym,
            Tn1.boot = pvalue.Tn1.boot,
            Tn2.boot = pvalue.Tn2.boot, x1 = x1, y1 = y1, x2 = x2, y2 = y2,
            m0x1.hat = m0x1.hat, m0x2.hat = m0x2.hat, eps1.hat = eps1.hat,
            eps01.hat = eps01.hat, eps2.hat = eps2.hat, eps02.hat = eps02.hat,
            m1x1.hat = m1x1.hat, m0x2.hat = m0x2.hat, sigma1x1.hat = sigma1x1.hat,
            sigma2x2.hat = sigma2x2.hat)

  class(r) <- c('list', 'function1')
  return(r)
}


#===============================================================================
# Example
#===============================================================================

n1 <- 100
n2 <- 200
x1 <- runif(n1)
y1 <- x1 + 0.25 * rnorm(n1)
x2 <- runif(n2)
y2 <- x2 + 0.25 * rnorm(n2)

res <- function1(x1, y1, x2, y2)





