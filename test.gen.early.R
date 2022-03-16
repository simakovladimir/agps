library(deSolve)

test.gen <- function(magnitude = list(rep(1, 3L), c(1, 1, 0.1)),
                     frequency = list(c(2 * pi, 2 * pi, 2 * pi), c(pi / 5, pi / 5, pi)),
                     phase     = list(rep(0, 3L), c(0, pi / 2, 0)),
                     sd        = rep(0, 19L),
                     q.0       = c(1, 0, 0, 0),
                     duration  = 10,
                     n         = 1L + 1E+3L,
                     g         = 9.80665,
                     E         = c(1.93254E+2, 2.4582E+1, 4.64090E+2)) {
  omega.fun <- as.function(alist(t = , magnitude[[1L]] * cos(drop(frequency[[1L]] %o% t) + phase[[1L]])))
  Omega.fun <- as.function(alist(q = , matrix(c(-q[2L], -q[3L], -q[4L],
                                                +q[1L], -q[4L], +q[3L],
                                                +q[4L], +q[1L], -q[2L],
                                                -q[3L], +q[2L], +q[1L]),
                                              4L, 3L, TRUE)))
  S.fun <- as.function(alist(q = ,
                             matrix(c(q[1L] ^ 2L + q[2L] ^ 2L - q[3L] ^ 2L - q[4L] ^ 2L,
                                      2 * (q[2L] * q[3L] + q[1L] * q[4L]),
                                      2 * (q[2L] * q[4L] - q[1L] * q[3L]),
                                      2 * (q[2L] * q[3L] - q[1L] * q[4L]),
                                      q[1L] ^ 2L - q[2L] ^ 2L + q[3L] ^ 2L - q[4L] ^ 2L,
                                      2 * (q[3L] * q[4L] + q[1L] * q[2L]),
                                      2 * (q[2L] * q[4L] + q[1L] * q[3L]),
                                      2 * (q[3L] * q[4L] - q[1L] * q[2L]),
                                      q[1L] ^ 2L - q[2L] ^ 2L - q[3L] ^ 2L + q[4L] ^ 2L),
                                    3L, 3L, TRUE)))
  q.dot.fun <- as.function(alist(t = , q = , p = , list(drop(0.5 * Omega.fun(q) %*% omega.fun(t)))))
  normalize <- as.function(alist(v = ,  v / sqrt(sum(v ^ 2L))))
  times <- seq(0, duration, length.out = n)
  omega <- omega.fun(times)
  q <- unname(cbind(q.0, sapply(2L : n - 1L, as.function(alist(k = ,
                                                               q.0 <<- normalize(lsoda(q.0,
                                                                                       times[0L : 1L + k],
                                                                                       q.dot.fun)[2L, -1L]))))))
  S <- lapply(1L : n, as.function(alist(k = , S.fun(q[, k]))))
  P <- +magnitude[[2L]] * frequency[[2L]] ^ 0L * cos(drop(frequency[[2L]] %o% times) + phase[[2L]])
  V <- -magnitude[[2L]] * frequency[[2L]] ^ 1L * sin(drop(frequency[[2L]] %o% times) + phase[[2L]])
  a <- -magnitude[[2L]] * frequency[[2L]] ^ 2L * cos(drop(frequency[[2L]] %o% times) + phase[[2L]]) +
    sapply(S, `%*%`, c(0, 0, g))
  B <- sapply(S, `%*%`, E)
  noise <- replicate(n, rnorm(19L, , sd))
  return(data.frame(time  = times,
                    P     = t(P + noise[1L : 3L, ]),
                    V     = t(V + noise[4L : 6L, ]),
                    q     = t(q + noise[7L : 10L, ]),
                    omega = t(omega + noise[11L : 13L, ]),
                    a     = t(a + noise[14L : 16L, ]),
                    B     = t(B + noise[17L : 19L, ])))
}

# library(rgl)
# plot3d(data$P.1, data$P.2, data$P.3, col = "navy", size = 3L, zlim = c(-1, 1), xlab = "", ylab = "", zlab = "")