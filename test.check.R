# WARNING!
# Temperature (to be included into control vector) is not implemented yet.

# State estimation kit for Assisted GPS navigation systems.
# Author: Vladimir Simakov (simakovladimir (at) gmail.com)
# Maintainer: Space Research Institute, Kiev, UA
# Dependencies: deSolve
# ...
# TODO:
# Enumerate and describe all functions (including arguments and algorithms):
#   cover such questions as how the arguments are organized, what rows and
#   cols mean in matrices and so on
# Describe all unobvious code episodes (e.g. serial update
#   in EKF implementation)
# Include into all functions appropriate validators (via stopifnot())
# Describe all wrappers for generic functions
# Describe specification of all envir.* functions including all variables
#   to be defined within them: x.dot.fun() -- process right side, h.fun()
#   -- observations right side,  normalize() -- either identity function,
#   or a function normalizing some part of the state vector (or entire
#   vector) during each iteration
# Note that NAs in measurements are supported
# Note that the name of the state variable within KF implementation should
#   always be `x`: this name is referred to in the output config utility to
#   display size of the state vectors set
# Signature of the KF methods should always be the same
# .rk4() -- one-step Runge-Kutta algorithm for KF predictor. If predictor is
#   not assumed (e.g. on parameters estimation) identity function may be used
#   instead. Moreover for prediction may be used such a construction as well:
#   normalize(unname(lsoda(x.0, time[k + 0L : 1L], x.dot.fun)[2L, -1L])),
#   where lsoda -- RK solver from deSolve package

library(deSolve)

.rk4 <- function(x, t, tau, fun) {
  
  k.1 <- fun(t, x)
  k.2 <- fun(t + tau / 2, x + tau * k.1 / 2)
  k.3 <- fun(t + tau / 2, x + tau * k.2 / 2)
  k.4 <- fun(t + tau, x + tau * k.3)
  
  return(x + tau * (k.1 + 2 * k.2 + 2 * k.3 + k.4) / 6)
  
# N.B. Eliminate dependency on t (see zhang2005.pdf) -- this is current
#      control input but not a function of t
#      Appropriate modifications within main code are required too
}

.aux.output.config <- function(env) {
  
  assign(
    "duration",
    proc.time()["elapsed"],
    envir = env)
  
  assign(
    "msg",
    c("complete\n",
      "A-GPS State Estimator\n",
      "Iteration %0*d/%0*d ",
      "Elapsed time: %d:%02d:%02d\n",
      "Grabbed size: %d bytes\n"),
    envir = env)
  
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
                    1L + !is.null(complete),
                    sprintf(env$msg[num], ...),
                    paste(
                      paste(
                        rep(".", 45L - nchar(sprintf(env$msg[num], ...))),
                        collapse = ""),
                      env$msg[complete])),
                  sep = "")),
              quote(flush.console())))))),
    envir = env)
  
  assign("oldfor", `for`,
         envir = env)
  
  assign(
    "for",
    switch(
      1L + env$output,
      env$oldfor,
      as.function(
        alist(
          loop.var = , loop.range = , loop.expr = ,
          do.call(
            `{`,
            list(
              N <- max(loop.range),
              fmtlen <- 1L + as.integer(trunc(log10(N))),
              env$msg.box(2L, NULL),
              eval(
                call(
                  "oldfor",
                  substitute(loop.var),
                  loop.range,
                  bquote(
                    do.call(
                      `{`,
                      list(
                        msg.box(3L, NULL, .(fmtlen), .(substitute(loop.var)),
                                .(fmtlen), .(N)),
                        .(substitute(loop.expr)),
                        msg.box(3L, 1L, .(fmtlen), .(substitute(loop.var)),
                                .(fmtlen), .(N)))))),
                envir = env),
              as.function(
                alist(
                  d = as.integer(round(proc.time()["elapsed"] - env$duration)),
                  env$msg.box(4L, NULL,
                              d %/% 3600L, (d %/% 60L) %% 60L, d %% 60L)))(),
              env$msg.box(5L, NULL, object.size(env$x))))))),
      envir = env)
}

ekf <- function(x.0, time, predict, ctrl, observ, sd, params,
                envir.deploy, output) {
  
  .aux.output.config(environment())
  envir.deploy(environment())
  
  n <- length(x.0)
  x <- x.0
  P <- diag(sd[[1L]] ^ 2L)
  Q <- diag(sd[[2L]] ^ 2L)
  R <- diag(sd[[3L]] ^ 2L)
  
  for (k in 2L : length(time) - 1L) {
    tau <- time[k + 1L] - time[k]
    x.0 <- normalize(unname(lsoda(x.0, time[k + 0L : 1L], x.dot.fun)[2L, -1L]))
    J <- J.fun(x.0, ctrl[k, ])
    G <- G.fun(x.0, ctrl[k, ])
    P <- `+`(
      (diag(n) + tau * J) %*% P %*% t(diag(n) + tau * J),
      tau ^ 2L * G %*% Q %*% t(G))
    H <- H.fun(x.0)
    y <- h.fun(x.0)
    z <- observ[k, ]
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
    x <- rbind(x, x.0, deparse.level = 0L)
  }
  
  return(x)
}

.agps.envir.deploy.generic <- function(env) {
  
  assign(
    "Omega.fun",
    as.function(
      alist(
        q = ,
        matrix(
          c(-q[2L], -q[3L], -q[4L],
            +q[1L], -q[4L], +q[3L],
            +q[4L], +q[1L], -q[2L],
            -q[3L], +q[2L], +q[1L]),
          4L, 3L, TRUE))),
    envir = env)
  
  assign(
    "S.fun",
    as.function(
      alist(
        q = ,
        matrix(
          c(q[1L] ^ 2L + q[2L] ^ 2L - q[3L] ^ 2L - q[4L] ^ 2L,
            2 * (q[2L] * q[3L] + q[1L] * q[4L]),
            2 * (q[2L] * q[4L] - q[1L] * q[3L]),
              2 * (q[2L] * q[3L] - q[1L] * q[4L]),
              q[1L] ^ 2L - q[2L] ^ 2L + q[3L] ^ 2L - q[4L] ^ 2L,
              2 * (q[3L] * q[4L] + q[1L] * q[2L]),
                2 * (q[2L] * q[4L] + q[1L] * q[3L]),
                2 * (q[3L] * q[4L] - q[1L] * q[2L]),
                q[1L] ^ 2L - q[2L] ^ 2L - q[3L] ^ 2L + q[4L] ^ 2L),
          3L, 3L, TRUE))),
    envir = env)
}

.agps.envir.deploy <- function(env) {
  
  .agps.envir.deploy.generic(env)
  
  assign(
    "u.k.fun",
    lapply(
      1L : 6L,
      as.function(
        alist(
          k = ,
          approxfun(env$time, env$ctrl[, k], rule = 2L)))),
    envir = env)
  
  assign(
    "u.fun",
    as.function(alist(t = , sapply(env$u.k.fun, do.call, list(t)))),
    envir = env)
  
  assign(
    "x.dot.fun",
    as.function(
      alist(
        t = , x = , p = ,
        list(
          c(
            x[4L : 6L],
            `+`(
              drop(t(env$S.fun(x[7L : 10L])) %*% env$u.fun(t)[4L : 6L]),
              c(0, 0, env$params$g)),
            drop(
              `%*%`(
                0.5 * env$Omega.fun(x[7L : 10L]),
                env$u.fun(t)[1L : 3L])))))),
    envir = env)
  
  assign(
    "h.fun",
    as.function(
      alist(
        x = ,
        c(x[1L : 6L],
          drop(env$S.fun(x[7L : 10L]) %*% env$params$E),
          `*`(
            env$params$p,
            exp(
              `/`(
                env$params$M * env$params$g * x[3L],
                env$params$R * env$params$t)))))),
    envir = env)
  
  assign(
    "normalize",
    as.function(
      alist(v = ,
            c(v[1L : 6L], v[7L : 10L] / sqrt(sum(v[7L : 10L] ^ 2L))))),
    envir = env)
  
  assign(
    "J.0.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (+x[07L] * u[4L] - x[10L] * u[5L] + x[09L] * u[6L]))),
    envir = env)
  
  assign(
    "J.1.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (+x[08L] * u[4L] + x[09L] * u[5L] + x[10L] * u[6L]))),
    envir = env)
  
  assign(
    "J.2.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (-x[09L] * u[4L] + x[08L] * u[5L] + x[07L] * u[6L]))),
    envir = env)
  
  assign(
    "J.3.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (-x[10L] * u[4L] - x[07L] * u[5L] + x[08L] * u[6L]))),
    envir = env)
  
  assign(
    "J.V.fun",
    as.function(
      alist(
        x = , u = ,
        matrix(
          c(
            +env$J.0.fun(x, u), +env$J.1.fun(x, u),
            +env$J.2.fun(x, u), +env$J.3.fun(x, u),
              -env$J.3.fun(x, u), -env$J.2.fun(x, u),
              +env$J.1.fun(x, u), +env$J.0.fun(x, u),
                +env$J.2.fun(x, u), -env$J.3.fun(x, u),
                -env$J.0.fun(x, u), +env$J.1.fun(x, u)),
          3L, 4L, TRUE))),
    envir = env)
  
  assign(
    "J.q.fun",
    as.function(
      alist(
        x = , u = ,
        0.5 * matrix(
          c(
            0.0E+0, -u[1L], -u[2L], -u[3L],
            +u[1L], 0.0E+0, +u[3L], -u[2L],
            +u[2L], -u[3L], 0.0E+0, +u[1L],
            +u[3L], +u[2L], -u[1L], 0.0E+0),
          4L, 4L, TRUE))),
    envir = env)
  
  assign(
    "J.fun",
    as.function(
      alist(
        x = , u = ,
        rbind(
          cbind(
            matrix(0, 3L, 3L),
            diag(3L),
            matrix(0, 3L, 4L)),
          cbind(
            matrix(0, 7L, 6L),
            rbind(env$J.V.fun(x, u), env$J.q.fun(x, u)))))),
    envir = env)
  
  assign(
    "G.fun",
    as.function(
      alist(
        x = , u = ,
        cbind(
          rbind(
            matrix(0, 6L, 3L),
            0.5 * env$Omega.fun(x[7L : 10L])),
          rbind(
            matrix(0, 3L, 3L),
            t(env$S.fun(x[7L : 10L])),
            matrix(0, 4L, 3L))))),
    envir = env)
  
  assign(
    "H.0.fun",
    as.function(
      alist(
        x = ,
        `*`(
          2,
          `+`(
            `+`(
              +x[07L] * env$params$E[1L],
              +x[10L] * env$params$E[2L]),
            -x[09L] * env$params$E[3L])))),
    envir = env)
  
  assign(
    "H.1.fun",
    as.function(
      alist(
        x = ,
        `*`(
          2,
          `+`(
            `+`(
              +x[08L] * env$params$E[1L],
              +x[09L] * env$params$E[2L]),
            +x[10L] * env$params$E[3L])))),
    envir = env)
  
  assign(
    "H.2.fun",
    as.function(
      alist(
        x = ,
        `*`(
          2,
          `+`(
            `+`(
              -x[09L] * env$params$E[1L],
              +x[08L] * env$params$E[2L]),
            -x[07L] * env$params$E[3L])))),
    envir = env)
  
  assign(
    "H.3.fun",
    as.function(
      alist(
        x = ,
        `*`(
          2,
          `+`(
            `+`(
              -x[10L] * env$params$E[1L],
              +x[07L] * env$params$E[2L]),
            +x[08L] * env$params$E[3L])))),
    envir = env)
  
  assign(
    "H.B.fun",
    as.function(
      alist(
        x = ,
        matrix(
          c(
            +env$H.0.fun(x), +env$H.1.fun(x),
            +env$H.2.fun(x), +env$H.3.fun(x),
              +env$H.3.fun(x), -env$H.2.fun(x),
              +env$H.1.fun(x), -env$H.0.fun(x),
                -env$H.2.fun(x), -env$H.3.fun(x),
                +env$H.0.fun(x), +env$H.1.fun(x)),
          3L, 4L, TRUE))),
    envir = env)
  
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
                `/`(
                  env$params$p * env$params$M * env$params$g,
                  env$params$R * env$params$t),
                exp(
                  `/`(
                    env$params$M * env$params$g * x[3L],
                    env$params$R * env$params$t))),
              rep(0, 3L))),
          rbind(
            matrix(0, 6L, 4L), env$H.B.fun(x), rep(0, 4L))))),
    envir = env)
}

agps.estimate <- function(data,
                          sd     = rep(0.1, 26L),
                          x.0    = as.numeric(data[1L, 2L : 11L]),
                          params = list(R = 8.3144598,
                                        M = 2.9E+4,
                                        g = 9.80665,
                                        E = c(1.93254E+2,
                                              2.45820E+1,
                                              4.64090E+2),
                                        p = 1.01325E+5,
                                        t = 2.9315E+2),
                          method = ekf,
                          output = TRUE) {
  
  x <- method(x.0,
              as.numeric(data$time),
              as.matrix(data[, 12L : 17L]),
              as.matrix(data[, c(2L : 7L, 18L : 21L)]),
              list(sd[01L : 10L],
                   sd[11L : 16L],
                   sd[17L : 26L]),
              params,
              .agps.envir.deploy,
              output)
  
  return(data.frame(time = data$time,
                    P    = x[, 1L : 3L],
                    V    = x[, 4L : 6L],
                    q    = x[, 7L : 10L]))
}

agps.simulate <- function(magnitude = list(rep(1, 3L),
                                           c(1, 1, 0.1)),
                          frequency = list(c(2 * pi, 2 * pi, 2 * pi),
                                           c(pi / 5, pi / 5, pi)),
                          phase     = list(rep(0, 3L),
                                           c(0, pi / 2, 0)),
                          sd        = rep(0, 20L),
                          q.0       = c(1, 0, 0, 0),
                          duration  = 10,
                          n         = 1L + 1E+3L,
                          params    = list(R = 8.3144598,
                                           M = 2.9E+4,
                                           g = 9.80665,
                                           E = c(1.93254E+2,
                                                 2.45820E+1,
                                                 4.64090E+2),
                                           p = 1.01325E+5,
                                           t = 2.9315E+2)) {
  
  .agps.envir.deploy.generic(environment())
  
  omega.fun <- as.function(
    alist(t = ,
          magnitude[[1L]] * cos(drop(frequency[[1L]] %o% t) + phase[[1L]])))

  q.dot.fun <- as.function(
    alist(t = , q = , p = ,
          list(drop(0.5 * Omega.fun(q) %*% omega.fun(t)))))
  
  normalize <- as.function(
    alist(v = ,
          v / sqrt(sum(v ^ 2L))))
  
  time <- seq(0, duration, length.out = n)
  omega <- omega.fun(time)
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
          c(0, 0, params$g)),
        2L, list),
      FALSE))
  B <- sapply(S, `%*%`, params$E)
  p <- params$p * exp(params$M * params$g * P[3L, ] / params$R / params$t)
  noise <- replicate(n, rnorm(20L, , sd))
  
  return(data.frame(time  = time,
                    P     = t(P + noise[1L : 3L, ]),
                    V     = t(V + noise[4L : 6L, ]),
                    q     = t(q + noise[7L : 10L, ]),
                    omega = t(omega + noise[11L : 13L, ]),
                    a     = t(a + noise[14L : 16L, ]),
                    B     = t(B + noise[17L : 19L, ]),
                    p     = p + noise[20L, ]))
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