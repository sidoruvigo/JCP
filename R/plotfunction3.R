
plot.function3(x, ...){
  par(mfrow = c(2, 2))

  plot(c(x1, x2),c(y1, y2), type = "n", xlab = "covariate", ylab = "response",
       main = "resgression functions")
  points(x1, y1, col = "red")
  points(x2, y2, col = "blue")
  lines(sort(x1), m1x1.hat[order(x1)], col = "red")
  lines(sort(x2), m2x2.hat[order(x2)], col = "blue")

  plot(c(x1, x2), c(y1, y2), type = "n", xlab = "covariate",
     ylab = "sigma", main = "conditional variance functions")
  lines(sort(x1), sigma1x1.hat[order(x1)], col = "red")
  lines(sort(x1), sigma0x1.hat[order(x1)], col = "red", lty = 2)
  lines(sort(x2), sigma2x2.hat[order(x2)], col = "blue")
  lines(sort(x2), sigma0x2.hat[order(x2)], col = "blue", lty = 2)

  plot(ecdf(eps1.hat), col = "red", main = "F_eps1")
  plot(ecdf(eps01.hat), col = "black", add = T)

  plot(ecdf(eps2.hat), col = "blue", main = "F_eps2")
  plot(ecdf(eps02.hat), col = "black", add = T)
}
