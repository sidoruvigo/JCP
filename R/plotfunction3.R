

plot.function3 <- function(x, ...) {
  graphics::par(mfrow = c(2, 2))

  graphics::plot(c(x$x1, x$x2), c(x$y1, x$y2), type = "n", xlab = "covariate", ylab = "response",
                 main = "regression functions")
  graphics::points(x$x1, x$y1, col = "red")
  graphics::points(x$x2, x$y2, col = "blue")
  graphics::lines(sort(x$x1), x$m1x1.hat[order(x$x1)], col = "red")
  graphics::lines(sort(x$x2), x$m2x2.hat[order(x$x2)], col = "blue")

  graphics::plot(c(x$x1, x$x2), c(x$y1, x$y2), type = "n", xlab = "covariate",
                 ylab = "sigma", main = "conditional variance functions")
  graphics::lines(sort(x$x1), x$sigma1x1.hat[order(x$x1)], col = "red")
  graphics::lines(sort(x$x1), x$sigma0x1.hat[order(x$x1)], col = "red", lty = 2)
  graphics::lines(sort(x$x2), x$sigma2x2.hat[order(x$x2)], col = "blue")
  graphics::lines(sort(x$x2), x$sigma0x2.hat[order(x$x2)], col = "blue", lty = 2)

  graphics::plot(ecdf(x$eps1.hat),  col = "red",   main = "F_eps1")
  graphics::plot(ecdf(x$eps01.hat), col = "black", add = T)

  graphics::plot(ecdf(x$eps2.hat),  col = "blue",  main = "F_eps2")
  graphics::plot(ecdf(x$eps02.hat), col = "black", add = T)
}
