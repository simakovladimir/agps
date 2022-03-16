# WARNING!
# Temperature (to be included into control vector) as well as
# pressure (to be included into observation vector) are not implemented yet.

library(deSolve)

state.gps <- function(data,
                      moments  = 1L : nrow(data),
                      sd       = rep(0.1, 19L),
                      x.0      = as.numeric(data[1L, 2L : 11L]),
                      duration = 10,
                      n        = 1L + 1E+3L,
                      g        = 9.80665,
                      E        = c(1.93254E+2, 2.4582E+1, 4.64090E+2)) {
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
  F.0.fun <- as.function(alist(x = , u = , 2 * (+x[07L] * u[4L] - x[10L] * u[5L] + x[09L] * u[6L])))
  F.1.fun <- as.function(alist(x = , u = , 2 * (+x[08L] * u[4L] + x[09L] * u[5L] + x[10L] * u[6L])))
  F.2.fun <- as.function(alist(x = , u = , 2 * (-x[09L] * u[4L] + x[08L] * u[5L] + x[07L] * u[6L])))
  F.3.fun <- as.function(alist(x = , u = , 2 * (-x[10L] * u[4L] - x[07L] * u[5L] + x[08L] * u[6L])))
  F.V.fun <- as.function(alist(x = , u = ,
                               matrix(c(+F.0.fun(x, u), +F.1.fun(x, u), +F.2.fun(x, u), +F.3.fun(x, u),
                                        -F.3.fun(x, u), -F.2.fun(x, u), +F.1.fun(x, u), +F.0.fun(x, u),
                                        +F.2.fun(x, u), -F.3.fun(x, u), -F.0.fun(x, u), +F.1.fun(x, u)),
                                      3L, 4L, TRUE)))
  F.q.fun <- as.function(alist(x = , u = ,
                               0.5 * matrix(c(     0, -u[1L], -u[2L], -u[3L],
                                              +u[1L],      0, +u[3L], -u[2L],
                                              +u[2L], -u[3L],      0, +u[1L],
                                              +u[3L], +u[2L], -u[1L],      0),
                                            4L, 4L, TRUE)))
  J.fun <- as.function(alist(x = , u = ,
                             rbind(cbind(matrix(0, 3L, 3L), diag(3L), matrix(0, 3L, 4L)),
                                   cbind(matrix(0, 7L, 6L), rbind(F.V.fun(x, u), F.q.fun(x, u))))))
  G.fun <- as.function(alist(x = , u = ,
                             cbind(rbind(matrix(0, 6L, 3L), 0.5 * Omega.fun(x[7L : 10L])),
                                   rbind(matrix(0, 3L, 3L), t(S.fun(x[7L : 10L])), matrix(0, 4L, 3L)))))
  H.0.fun <- as.function(alist(x = , y = , 2 * (+x[07L] * y[7L] + x[10L] * y[8L] - x[09L] * y[9L])))
  H.1.fun <- as.function(alist(x = , y = , 2 * (+x[08L] * y[7L] + x[09L] * y[8L] + x[10L] * y[9L])))
  H.2.fun <- as.function(alist(x = , y = , 2 * (-x[09L] * y[7L] + x[08L] * y[8L] - x[07L] * y[9L])))
  H.3.fun <- as.function(alist(x = , y = , 2 * (-x[10L] * y[7L] + x[07L] * y[8L] + x[08L] * y[9L])))
  H.B.fun <- as.function(alist(x = , y = ,
                               matrix(c(+H.0.fun(x, y), +H.1.fun(x, y), +H.2.fun(x, y), +H.3.fun(x, y),
                                        +H.3.fun(x, y), -H.2.fun(x, y), +H.1.fun(x, y), -H.0.fun(x, y),
                                        -H.2.fun(x, y), -H.3.fun(x, y), +H.0.fun(x, y), +H.1.fun(x, y)),
                                      3L, 4L, TRUE)))
  H.fun <- as.function(alist(x = , y = ,
                             cbind(rbind(diag(6L), matrix(0, 3L, 6L)),
                                   rbind(matrix(0, 6L, 4L), cbind(H.B.fun(x, y))))))
  u.k.fun <- lapply(12L : 17L, as.function(alist(k = , approxfun(data$time, data[, k], rule = 2L))))
  u.fun <- as.function(alist(t = , sapply(u.k.fun, do.call, list(t))))
  x.dot.fun <- as.function(alist(t = , x = , p = ,
                                 list(c(x[4L : 6L],
                                        drop(t(S.fun(x[7L : 10L])) %*% u.fun(t)[4L : 6L]) - c(0, 0, g),
                                        drop(0.5 * Omega.fun(x[7L : 10L]) %*% u.fun(t)[1L : 3L])))))
  normalize <- as.function(alist(v = , c(v[1L : 6L], v[7L : 10L] / sqrt(sum(v[7L : 10L] ^ 2L)))))
  P <- diag(sd[1L : 10L] ^ 2L)
  Q <- diag(sd[11L : 16L] ^ 2L)
  R <- diag(sd[c(1L : 6L, 17L : 19L)] ^ 2L)
  tau <- duration / (n - 1L)
  times <- seq(0, duration, length.out = n)
  x <- x.0
  for (k in moments[-1L]) {
    x.0 <- normalize(unname(lsoda(x.0, times[0L : 1L + k], x.dot.fun)[2L, -1L]))
    J <- J.fun(x.0, as.numeric(data[k, 12L : 17L]))
    G <- G.fun(x.0, as.numeric(data[k, 12L : 17L]))
    P <- (diag(10L) + tau * J) %*% P %*% t(diag(10L) + tau * J) + tau ^ 2L * G %*% Q %*% t(G)
    H <- H.fun(x.0, as.numeric(data[k, c(2L : 7L, 18L : 20L)]))
    y <- c(x.0[1L : 6L], drop(S.fun(x.0[7L : 10L]) %*% E))
    K <- P %*% t(H) %*% chol2inv(chol(H %*% P %*% t(H) + R))
    x.0 <- normalize(x.0 + drop(K %*% (as.numeric(data[k, c(2L : 7L, 18L : 20L)]) - y)))
    x <- rbind(x, x.0, deparse.level = 0L)
    P <- P - K %*% H %*% t(P)
    cat(k, ": ok\n")
    flush.console()
  }
  return(data.frame(time = times[moments],
                    P    = x[, 1L : 3L],
                    V    = x[, 4L : 6L],
                    q    = x[, 7L : 10L]))
}