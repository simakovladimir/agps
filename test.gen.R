library(deSolve)

test.gen <- function(magnitude = list(rep(1, 3L),
                                      c(1, 1, 0.1)),
                     frequency = list(c(2 * pi, 2 * pi, 2 * pi),
                                      c(pi / 5, pi / 5, pi)),
                     phase     = list(rep(0, 3L),
                                      c(0, pi / 2, 0)),
                     sd        = rep(0, 20L),
                     q.0       = c(1, 0, 0, 0),
                     duration  = 10,
                     n         = 1L + 1E+3L,
                     env       = list(R = 8.3144598,
                                      M = 2.9E+4,
                                      g = 9.80665,
                                      E = c(1.93254E+2, 2.45820E+1, 4.64090E+2),
                                      p = 1.01325E+5,
                                      t = 2.9315E+2)) {
  assign(
    "omega.fun",
    as.function(
      alist(
        t = ,
        magnitude[[1L]] * cos(drop(frequency[[1L]] %o% t) + phase[[1L]]))))
  assign(
    "Omega.fun",
    as.function(
      alist(
        q = ,
        matrix(
          c(
            -q[2L], -q[3L], -q[4L],
            +q[1L], -q[4L], +q[3L],
            +q[4L], +q[1L], -q[2L],
            -q[3L], +q[2L], +q[1L]),
          4L, 3L, TRUE))))
  assign(
    "S.fun",
    as.function(
      alist(
        q = ,
        matrix(
          c(
            q[1L] ^ 2L + q[2L] ^ 2L - q[3L] ^ 2L - q[4L] ^ 2L,
            2 * (q[2L] * q[3L] + q[1L] * q[4L]),
            2 * (q[2L] * q[4L] - q[1L] * q[3L]),
              2 * (q[2L] * q[3L] - q[1L] * q[4L]),
              q[1L] ^ 2L - q[2L] ^ 2L + q[3L] ^ 2L - q[4L] ^ 2L,
              2 * (q[3L] * q[4L] + q[1L] * q[2L]),
                2 * (q[2L] * q[4L] + q[1L] * q[3L]),
                2 * (q[3L] * q[4L] - q[1L] * q[2L]),
                q[1L] ^ 2L - q[2L] ^ 2L - q[3L] ^ 2L + q[4L] ^ 2L),
          3L, 3L, TRUE))))
  assign(
    "q.dot.fun",
    as.function(
      alist(
        t = , q = , p = ,
        list(
          drop(0.5 * Omega.fun(q) %*% omega.fun(t))))))
  assign(
    "normalize",
    as.function(
      alist(v = ,
            v / sqrt(sum(v ^ 2L)))))
  
  time <- seq(0, duration, length.out = n)
  omega <- omega.fun(time)
  
#   q <- t(normalize(unname(lsoda(q.0, times, q.dot.fun)[, -1L])))
  
  q <- unname(
    cbind(
      q.0,
      sapply(
        2L : n - 1L,
        as.function(
          alist(
            k = ,
            normalize(
              lsoda(
                q.0,
                time[0L : 1L + k],
                q.dot.fun)[2L, -1L]) ->> q.0)))))
  S <- lapply(1L : n, as.function(alist(k = , S.fun(q[, k]))))
  P <- `*`(magnitude[[2L]],
           cos(drop(frequency[[2L]] %o% time) + phase[[2L]]))
  V <- `*`(-magnitude[[2L]] * frequency[[2L]],
           sin(drop(frequency[[2L]] %o% time) + phase[[2L]]))
  a <- mapply(
    `%*%`, S,
    unlist(
      apply(
        `-`(
          `*`(
            -magnitude[[2L]] * frequency[[2L]] ^ 2L,
            cos(drop(frequency[[2L]] %o% time) + phase[[2L]])),
          c(0, 0, env$g)),
        2L, list),
      FALSE))
  B <- sapply(S, `%*%`, env$E)
  p <- env$p * exp(env$M * env$g * P[3L, ] / env$R / env$t)
  noise <- replicate(n, rnorm(20L, , sd))
  
# #   return(c(approxfun(times, omega[1L, ], rule = 2L)(0), approxfun(times, omega[1L, ], rule = 2L)(times[n - 1L]), approxfun(times, omega[1L, ], rule = 2L)(times[n]), approxfun(times, omega[1L, ], rule = 2L)(11)))
#   
#   q.0 <- c(1, 0, 0, 0)
#   q.est <- unname(
#     cbind(
#       q.0,
#       sapply(
#         2L : n - 1L,
#         as.function(
#           alist(
#             k = ,
#             normalize(
#               lsoda(
#                 q.0,
#                 times[0L : 1L + k],
#                 as.function(
#                   alist(
#                     t = , q = , p = ,
#                     list(
#                       drop(
#                         `%*%`(
#                           0.5 * Omega.fun(q),
#                           apply(
#                             omega, 1L,
#                             as.function(
#                               alist(
#                                 v = ,
#                                 approxfun(times, v, rule = 2L)(t))))))))))[2L, -1L]) ->> q.0)))))
# #   return(list(
# #     original = q,
# #     estimated = q.est))
#   S.est <- lapply(1L : n, as.function(alist(k = , t(S.fun(q.est[, k])))))
#   a.est <- mapply(`%*%`, S.est, unlist(apply(a, 2L, list), FALSE)) - c(0, 0, g)
# #   return(list(
# #     original = -magnitude[[2L]] * frequency[[2L]] ^ 2L * cos(drop(frequency[[2L]] %o% times) + phase[[2L]]),
# #     estimated = a.est))
  
  return(data.frame(time  = time,
                    P     = t(P + noise[1L : 3L, ]),
                    V     = t(V + noise[4L : 6L, ]),
                    q     = t(q + noise[7L : 10L, ]),
                    omega = t(omega + noise[11L : 13L, ]),
                    a     = t(a + noise[14L : 16L, ]),
                    B     = t(B + noise[17L : 19L, ]),
                    p     = p + noise[20L, ]))
}

# library(rgl)
# plot3d(data$P.1, data$P.2, data$P.3, col = "navy", size = 3L, zlim = c(-1, 1), xlab = "", ylab = "", zlab = "")
# plot3d(x$P.1[1L : 500L], x$P.2[1L : 500L], x$P.3[1L : 500L], col = "red", size = 5L, add = TRUE)
# plot3d(x$P.1[500L : 750L], x$P.2[500L : 750L], x$P.3[500L : 750L], col = "green", size = 5L, add = TRUE)