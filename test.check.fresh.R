# WARNING!
# Temperature (to be included into control vector) is not implemented yet.

library(deSolve)

state.gps <- function(data,
                      sd  = rep(0.1, 20L),
                      x.0 = as.numeric(data[1L, 2L : 11L]),
                      env = list(R = 8.3144598,
                                 M = 2.9E+4,
                                 g = 9.80665,
                                 E = c(1.93254E+2, 2.4582E+1, 4.64090E+2),
                                 p = 1.01325E+5,
                                 t = 2.9315E+2)) {
  
  duration <- proc.time()["elapsed"]
  
  msg<-c(
    "Initialization ",
    "Iteration %0*d/%0*d ",
    "Elapsed time: %d:%02d:%02d\n",
    "Grabbed size: %d bytes\n")
  assign(
    "msg.box",
    as.function(
      alist(
        num = , complete = , ... = ,
        invisible(
          do.call(
            `{`,
            list(
              quote(
                cat(
                  switch(
                    1L + complete,
                    sprintf(msg[num], ...),
                    paste(
                      paste(
                        rep(".", 45L - nchar(sprintf(msg[num], ...))),
                        collapse = ""),
                      "complete\n")),
                  sep = "")),
              quote(flush.console())))))))
  
  msg.box(1L, FALSE)
  
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
    "F.0.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (+x[07L] * u[4L] - x[10L] * u[5L] + x[09L] * u[6L]))))
  assign(
    "F.1.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (+x[08L] * u[4L] + x[09L] * u[5L] + x[10L] * u[6L]))))
  assign(
    "F.2.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (-x[09L] * u[4L] + x[08L] * u[5L] + x[07L] * u[6L]))))
  assign(
    "F.3.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (-x[10L] * u[4L] - x[07L] * u[5L] + x[08L] * u[6L]))))
  assign(
    "F.V.fun",
    as.function(
      alist(
        x = , u = ,
        matrix(
          c(
            +F.0.fun(x, u), +F.1.fun(x, u),
            +F.2.fun(x, u), +F.3.fun(x, u),
              -F.3.fun(x, u), -F.2.fun(x, u),
              +F.1.fun(x, u), +F.0.fun(x, u),
                +F.2.fun(x, u), -F.3.fun(x, u),
                -F.0.fun(x, u), +F.1.fun(x, u)),
          3L, 4L, TRUE))))
  assign(
    "F.q.fun",
    as.function(
      alist(
        x = , u = , 0.5 * matrix(c(0.0E+0, -u[1L], -u[2L], -u[3L],
                                   +u[1L], 0.0E+0, +u[3L], -u[2L],
                                   +u[2L], -u[3L], 0.0E+0, +u[1L],
                                   +u[3L], +u[2L], -u[1L], 0.0E+0),
                                 4L, 4L, TRUE))))
  assign(
    "J.fun",
    as.function(
      alist(
        x = , u = ,
        rbind(
          cbind(
            matrix(0, 3L, 3L), diag(3L), matrix(0, 3L, 4L)),
          cbind(
            matrix(0, 7L, 6L), rbind(F.V.fun(x, u), F.q.fun(x, u)))))))
  assign(
    "G.fun",
    as.function(
      alist(
        x = , u = ,
        cbind(
          rbind(
            matrix(0, 6L, 3L), 0.5 * Omega.fun(x[7L : 10L])),
          rbind(
            matrix(0, 3L, 3L), t(S.fun(x[7L : 10L])), matrix(0, 4L, 3L))))))
  assign(
    "H.0.fun",
    as.function(
      alist(
        x = ,
        2 * (+x[07L] * env$E[1L] + x[10L] * env$E[2L] - x[09L] * env$E[3L]))))
  assign(
    "H.1.fun",
    as.function(
      alist(
        x = ,
        2 * (+x[08L] * env$E[1L] + x[09L] * env$E[2L] + x[10L] * env$E[3L]))))
  assign(
    "H.2.fun",
    as.function(
      alist(
        x = ,
        2 * (-x[09L] * env$E[1L] + x[08L] * env$E[2L] - x[07L] * env$E[3L]))))
  assign(
    "H.3.fun",
    as.function(
      alist(
        x = ,
        2 * (-x[10L] * env$E[1L] + x[07L] * env$E[2L] + x[08L] * env$E[3L]))))
  assign(
    "H.B.fun",
    as.function(
      alist(
        x = ,
        matrix(
          c(
            +H.0.fun(x), +H.1.fun(x), +H.2.fun(x), +H.3.fun(x),
            +H.3.fun(x), -H.2.fun(x), +H.1.fun(x), -H.0.fun(x),
            -H.2.fun(x), -H.3.fun(x), +H.0.fun(x), +H.1.fun(x)),
          3L, 4L, TRUE))))
  assign(
    "H.fun",
    as.function(
      alist(
        x = ,
        cbind(
          rbind(
            diag(6L),
            matrix(0, 3L, 6L),
            c(
              rep(0, 2L),
              `*`(
                env$p * env$M * env$g / env$R / env$t,
                exp(env$M * env$g * x[3L] / env$R / env$t)),
              rep(0, 3L))),
          rbind(
            matrix(0, 6L, 4L), H.B.fun(x), rep(0, 4L))))))
#   assign(
#     "H.fun",
#     as.function(
#       alist(
#         x = ,
#         cbind(
#           rbind(
#             matrix(0, 3L, 6L),
#             c(
#               rep(0, 2L),
#               `*`(
#                 env$p * env$M * env$g / env$R / env$t,
#                 exp(env$M * env$g * x[3L] / env$R / env$t)),
#               rep(0, 3L))),
#           rbind(
#             H.B.fun(x), rep(0, 4L))))))
#   assign(
#     "H.fun",
#     as.function(
#       alist(
#         x = ,
#         cbind(matrix(0, 3L, 6L), H.B.fun(x)))))
  assign(
    "u.k.fun",
    lapply(
      12L : 17L,
      as.function(
        alist(
          k = ,
          approxfun(data$time, data[, k], rule = 2L)))))
  assign(
    "u.fun",
    as.function(alist(t = , sapply(u.k.fun, do.call, list(t)))))
  assign(
    "x.dot.fun",
    as.function(
      alist(
        t = , x = , p = ,
        list(
          c(
            x[4L : 6L],
            `+`(
              drop(t(S.fun(x[7L : 10L])) %*% u.fun(t)[4L : 6L]),
              c(0, 0, env$g)),
            drop(0.5 * Omega.fun(x[7L : 10L]) %*% u.fun(t)[1L : 3L]))))))
  assign(
    "normalize",
    as.function(
      alist(v = ,
            c(v[1L : 6L], v[7L : 10L] / sqrt(sum(v[7L : 10L] ^ 2L))))))
  
  P <- diag(sd[1L : 10L] ^ 2L)
  Q <- diag(sd[11L : 16L] ^ 2L)
  R <- diag(sd[c(1L : 6L, 17L : 20L)] ^ 2L)
#   R <- diag(sd[c(17L : 20L)] ^ 2L)
#   R <- diag(sd[c(17L : 19L)] ^ 2L)
  n <- nrow(data) - 1L
  fmtlen <- as.integer(trunc(log10(n))) + 1L
  x <- x.0
  
  msg.box(1L, TRUE)
  
  for (k in 1L : n) {
    msg.box(2L, FALSE, fmtlen, k, fmtlen, n)
    tau <- data$time[k + 1L] - data$time[k]
    x.0 <- normalize(unname(lsoda(x.0, data$time[0L : 1L + k],
                                  x.dot.fun)[2L, -1L]))
    J <- J.fun(x.0, as.numeric(data[k, 12L : 17L]))
    G <- G.fun(x.0, as.numeric(data[k, 12L : 17L]))
    P <- `+`(
      (diag(10L) + tau * J) %*% P %*% t(diag(10L) + tau * J),
      tau ^ 2L * G %*% Q %*% t(G))
    H <- H.fun(x.0)
#     y <- drop(H %*% x.0)
    y <- c(x.0[1L : 6L],
           drop(S.fun(x.0[7L : 10L]) %*% env$E),
           env$p * exp(env$M * env$g * x.0[3L] / env$R / env$t))
#     y <- c(drop(S.fun(x.0[7L : 10L]) %*% env$E),
#            env$p * exp(env$M * env$g * x.0[3L] / env$R / env$t))
#     y <- drop(S.fun(x.0[7L : 10L]) %*% env$E)
    z <- as.numeric(data[k, c(2L : 7L, 18L : 21L)])
#     z <- as.numeric(data[k, c(18L : 21L)])
#     z <- as.numeric(data[k, c(18L : 20L)])
    R.k <- R[!is.na(z), !is.na(z)]
    H.k <- H[!is.na(z), ]
    y.k <- y[!is.na(z)]
    z.k <- z[!is.na(z)]
    x.0 <- normalize(
      `+`(
        x.0,
        as.function(
          alist(
            i = 1L : nrow(R.k), U = P,
            switch(
              1L + (length(i) != 0),
              0 * U[1L, 1L],
              `+`(
                `*`(
                  `/`(
                    drop(U %*% H.k[i[1L], ]),
                    `+`(
                      drop(H.k[i[1L], ] %*% U %*% H.k[i[1L], ]),
                      R.k[i[1L], i[1L]])) -> K,
                  z.k[i[1L]] - y.k[i[1L]]),
                sys.function(0L)(
                  i[-1L],
                  U - as.matrix(K) %*% H.k[i[1L], ] %*% U ->> P)))))()))
#     x.0 <- normalize(
#       `+`(
#         x.0,
#         as.function(
#           alist(
#             i = 1L : nrow(R), U = P, x = x.0,
#             switch(
#               1L + (length(i) != 0),
#               0 * U[1L, 1L],
#               `+`(
#                 `*`(
#                   `/`(
#                     drop(U %*% H[i[1L], ]),
#                     `+`(
#                       drop(H[i[1L], ] %*% U %*% H[i[1L], ]),
#                       R[i[1L], i[1L]])) -> K,
#                   z[i[1L]] - drop(H[i[1L], ] %*% x)) -> add,
#                 sys.function(0L)(
#                   i[-1L],
#                   U - as.matrix(K) %*% H[i[1L], ] %*% U ->> P,
#                   x + add)))))()))
    x <- rbind(x, x.0, deparse.level = 0L)
    msg.box(2L, TRUE, fmtlen, k, fmtlen, n)
  }
  
  duration <- as.integer(round(proc.time()["elapsed"] - duration))
  msg.box(3L, FALSE, duration%/%3600L, (duration %/% 60L) %% 60L, duration %% 60L)
  msg.box(4L, FALSE, object.size(x))
  
  return(data.frame(time = data$time,
                    P    = x[, 1L : 3L],
                    V    = x[, 4L : 6L],
                    q    = x[, 7L : 10L]))
}

# N.B. Accuracy of quaternion in different cases:
# Barometer, GPS
# > range(x$q.1 - data$q.1)
# [1] -0.004336159  0.002784675
# > range(x$q.2 - data$q.2)
# [1] -0.005881168  0.014526472
# > range(x$q.3 - data$q.3)
# [1] -0.005098299  0.006327247
# > range(x$q.4 - data$q.4)
# [1] -0.007814326  0.033279578
# 
# Barometer, No GPS
# > range(x$q.1 - data$q.1)
# [1] -0.007985431  0.004317189
# > range(x$q.2 - data$q.2)
# [1] -0.02738764  0.01065087
# > range(x$q.3 - data$q.3)
# [1] -0.007569877  0.005766134
# > range(x$q.4 - data$q.4)
# [1] -0.05611588  0.01904716
# 
# No Barometer, No GPS
# > range(x$q.1 - data$q.1)
# [1] -0.009582076  0.004805754
# > range(x$q.2 - data$q.2)
# [1] -0.03145091  0.01090315
# > range(x$q.3 - data$q.3)
# [1] -0.008171277  0.005792555
# > range(x$q.4 - data$q.4)
# [1] -0.06479341  0.01909136