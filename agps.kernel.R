# Copyright 2015 Vladimir Simakov
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# AGPS Intermediate Notes
# LAST UPDATE: 2015/11/03

# WARNING!
# Temperature (to be included into control vector) is not implemented yet.

# BRIEF SUMMARY

# State estimation kit for Assisted GPS (AGPS) navigation systems
# Author: Vladimir Simakov
# Dependencies: r-base
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
#   instead. Moreover for prediction such a construction may be used as well:
#   normalize(unname(lsoda(x.0, time[k + 0L : 1L], x.dot.fun)[2L, -1L])),
#   where lsoda -- RK solver from deSolve package (fragment above with lsoda
#   is incorrect, validate it before publishing)
# 
# N.B. For barometer altitude measurements it is possible to implement
#   measurement expression as log(p + w), instead of pure p + w, where w --
#   white noise. It allows to simplify calculations and avoid huge numbers
#   (close to infinity) that damage calculations. Further -- to linearize by
#   w and obtain "almost" identity matrix.
# 
# N.B. Implement agps.simulate() without referring to numerical procedures,
#   i.e. purely analytically. It would allow to instantly obtain normalized
#   solution since due to properties of the equation all its solutions conserve
#   norm of the initial vector. UPDATE: Analytical Solver is too hard to
#   implement as the equation for q has periodic coefficients. There is special
#   theory which allows to build general solution for such systems, however
#   numerical solutions are absolutely applicable here. (See Floquet Theory,
#   in particular: R. Grinshaw, Nonlinear ODEs.)
#
# Accurately specify references ([]-format) within comments here.
#
# For .cholupdate() put reference to arXiv article on low rank Cholesky update.
#   Note: .cholupdate() replaced with .schud().
#
# N.B. All operations of SR-UKF are well implemented within this code. The only
#   serious issue referes to the Cholesky one-rank downdate. MATLAB function
#   cholupdate() discussed in Wiki and implemented here is dramatically
#   unstable and needs to be modified. See Dongarra et al, p. 182 to learn more
#   about algorithmic details of such a modification.
#
# .schud() and .schdd() are appropriate functions for Cholesky Decomposition
#   rank-one update and downdate (respectively). Here is an implementation of
#   simplified and numerically stable algorithms originated from LINPACK.
#   Note that .schud() assumes the last argument to be strictly positive.
#   Although it is allowed to be negative (downdating mode), this case leads
#   to stability losses. It is recommended to use .schdd() instead.
#   Note also that both functions are designed to work with lower triangular
#   matrices (in contrast to LINPACK versions where upper triangulars are 
#   used). It is also expected that these functions will be applied for vectors
#   and matrices whose dimension is greater than 1. All necessary validity
#   tests are omitted, user should check data before applying the functions.
#
# Note that there are two equivalent implementations of .schud() and .schdd()
#   functions. The first pair (.schud.rec() and .schdd.rec()) is implemented
#   via recursion. The second pair (.schud() and .schdd()) is implemented via
#   functions of apply family. Both pairs give identical results and have the
#   same performance. Although specified functions give indistinguishible
#   results, it is recommended to prefer apply-based ones. The recurrent
#   implementations may be used as an example of recursive constructions in R
#   (however they also may be replaced with basic Reduce() function).
#
# Note also that within "method" functions (ekf(), ukf(), etc.) .schud() and
#   .schdd() functions should be referred to by means of names cholupdate()
#   and choldowndate() respectively. Appropriate assignments should be
#   performed within the deployment routines.

# In ukf.sr(), to compare results of square-root modification with standard
#   implementation of UKF, run at the end of main cycle (_before_ the last
#   S[1L : n, 1L : n] assignment) the followng code:
#     S.x <- t(chol(`-`(S[1L : n, 1L : n] %*% t(S[1L : n, 1L : n])),
#                   K %*% (S.y %*% t(S.y)) %*% t(K)))
#   and compare results given by Reduce(.schdd, ...) function for
#   S[1L : n, 1L : n] with the result placed in S.x. They should be close to
#   each other.
#
# According to initial experiments, the difference in quality between UKF and
#   EKF is almost absent. They both give the results on the similar level of
#   accuracy even on partially empty data sets. Moreover, EKF is being executed
#   much faster than UKF.
#
# TODO: Unify implementation of ekf() and ukf() making it independent on
#   problem definition (exclude references to A-GPS variables). Make unified
#   subroutines for the right side and so on.
#   Also think about problems with output when .aux.output.config() is placed
#   out of method() function.
#
# UPDATE: normalize() is replaced by adjust().
#
# Within the method() functions (ekf(), ukf()) only one (!) iteration should be
# implemented. Cycles are not presumed.
#
# Since process covariance is to be updated on each iteration as well as the
#   state vector, method functions (ekf(), ukf()) are designed to take as an
#   argument and to return on exit a list contatining both the renewed x and
#   P (state vector and covariance matrix respectively).
#
# The very first argument of the method functions (ekf(), ukf()) is compounded
#   of two data structures: state vector and information that allows to
#   reproduce current covariance matrix. The latter is either vector of
#   variances (components of the state vectors are assumed to be independent),
#   or current estimation of covariance matrix, or its square root (in the
#   Cholesky meaning, for algorithms referring to UKF family). All these cases
#   should be considered. For performance reasons method functions don't carry
#   out any validity checks. Initial value of argument xP should consist of
#   a priory estimate of the state vector and vector with variances for it.
#   Such a format allows method function to distinct first run from subsequent
#   runs within external loop (implemented via the Reduce() function).
#
# Return value of method functions must be compatible with the first argument.
#
# UPDATE. In schud(), removed ability to update decomposition with the scaled
#   vector. It does not have actual sense. Appropriate code fragment has been
#   modified.
#
# N.B. Having estimated parameters of magnetometer, put them into the returned
#   result of the function agps.estimate().
#
# TODO:
# Cast the following Reduce() results to matrix form (not list of list as it
#   is right now). Think on how to add into the resulting matrix estimations
#   from magnetometer calibration: either additional element of entire list,
#   or additional three components of state vector (which of these whould be
#   better and faster and makes code clearer). Possible solution for current
#   problem definition: sapply(Reduce(...), `[[`, 1L). This should give the
#   expected behavior, but in the absence of magnetometer corrections. When
#   magnetometer parameters are added, this works only if we include them as
#   additional components straight to the state vector. [DONE]
#
# TODO: In the current rk4() implementaton the three-dot argument placed
#   inside a subfunction does not properly work. Think on how to fix.
#   [FIXED: ACTUALLY, IT'S NOT AN ISSUE]
#
# TODO: Loop-less implementation works properly, time losses are not observed
#   (that is especially pleasant). However, ekf() + Reduce() gives strange
#   result. Check it in detail.
#
# CHALLENGE COMPLETE. Loop-less implementation is debugged and works fine.
#
# TODO. Implement parameter estimating for magnetometer calibration using
#   existing method() functions. See Van der Merwe paper for details. Use
#   Soft Iron / Hard Iron model for appropriate measurement function (AN4246).
#
# TODO.
#   (i) Specify Greenwich as the default magnetic region. [DONE]
#   (ii) Estimate magnetic parameters as with the use of currently obtained
#     quaternion orientation, as by means of autonomous norm fitting (see
#     the remark above).
#   (iii) Make a use of bias models for both gyro and accelerometer values.
#   (iv) Include additional error variances for recently considered errors.
#   (v) Introduce formula for altitude measurement via the pressure values.
#
# N.B. Important observation on the R behavior.
#   When using functional-style subset operator `[`, it is parsed in specific
#   manner. E.g. `[`(a + b, i) is being expanded into a + b[i]. This means that
#   subset operator is applied to the unparsed string of an expression, but not
#   to the parsed result of its evaluation. It should be considered in order to
#   avoid unexpected logical errors. Thus, for recently studied case one should
#   surround "internal" expression with brackets `(` to obtain expected result.
#   `[`((a + b), i) is the same as `[`(`(`(`+`(a, b)), i)) and gives (a + b)[i]

ekf <- function(x, P, Q, R, u, z, x.fun, z.fun,
                ekf.J.fun, ekf.G.fun, ekf.H.fun, ...) {
  
  P <- switch(1L + is.null(dim(P)), P, diag(P ^ 2L))
  Q <- diag(Q ^ 2L)
  R <- diag(R ^ 2L)
  
  n <- length(x)
  m <- length(u)
  s <- length(z)
  
  x <- x.fun(x, u, numeric(m))
  J <- ekf.J.fun(x, u)
  G <- ekf.G.fun(x, u)
  P <- (diag(n) + J) %*% P %*% t(diag(n) + J) + G %*% Q %*% t(G)
  H <- ekf.H.fun(x)
  y <- z.fun(x, numeric(s), u)
  xP <- Reduce(
    as.function(
      alist(
        v = , i = ,
        list(
          `+`(`*`(`/`(drop(v[[2L]] %*% H[i, ]),
                      `+`(drop(H[i, ] %*% v[[2L]] %*% H[i, ]),
                          R[i, i])) -> K,
                  z[i] - y[i]),
              v[[1L]]),
          v[[2L]] - as.matrix(K) %*% H[i, ] %*% v[[2L]]))),
    seq_along(z), list(x, P))
  
  return(xP)
}

ukf.sr <- function(x, P, Q, R, u, z, x.fun, z.fun,
                   ukf.alpha = 1, ukf.beta = 2, ukf.kappa = 0,
                   ukf.chud, ukf.chdd, ...) {
  
  P <- switch(1L + is.null(dim(P)), P, diag(P))
  Q <- diag(Q)
  R <- diag(R)
  
  n <- length(x)
  m <- length(u)
  s <- length(z)
  L <- n + m + s
  
  lambda <- ukf.alpha ^ 2L * (L + ukf.kappa) - L
  gamma <- sqrt(L + lambda)
  w.m <- lambda / (L + lambda)
  w.c <- w.m + 1 - ukf.alpha ^ 2L + ukf.beta
  w <- 1 / 2 / (L + lambda)
  X <- unname(
    `+`(
      c(x, numeric(m + s)),
      Reduce(
        cbind,
        lapply(
          c(1, -1), `*`,
          `*`(
            rbind(
              cbind(P, matrix(0, n, m + s)),
              cbind(matrix(0, m, n), Q, matrix(0, m, s)),
              cbind(matrix(0, s, m + n), R)),
            gamma)), numeric(L))))
  X[1L : n, ] <- apply(
    X[1L : (n + m), ], 2L,
    as.function(
      alist(
        v = ,
        x.fun(v[1L : n], u, v[n + 1L : m]))))
  x <- w.m * X[1L : n, 1L] + w * apply(X[1L : n, -1L], 1L, sum)
  P <- ukf.chud(
    t(qr.R(qr(sqrt(w) * t(X[1L : n, -1L] - x)))),
    sqrt(w.c) * (X[1L : n, 1L] - x))
  xP <- switch(
    1L + (length(z) != 0L), list(x, P),
    `{`(
      Y <- mapply(z.fun,
                  lapply(apply(X[1L : n, ], 2L, list), unlist),
                  lapply(apply(X[n + m + 1L : s, ], 2L, list), unlist),
                  MoreArgs = list(u)),
      y <- w.m * Y[, 1L] + w * apply(Y[, -1L], 1L, sum),
      P.y <- ukf.chud(
        t(qr.R(qr(sqrt(w) * t(Y[, -1L] - y)))),
        sqrt(w.c) * (Y[, 1L] - y)),
      P.x.y <- `+`(
        w.c * ((X[1L : n, 1L] - x) %o% (Y[, 1L] - y)),
        `*`(
          Reduce(
            `+`,
            Map(
              `%o%`,
              lapply(apply(X[1L : n, -1L] - x, 2L, list), unlist),
              lapply(apply(Y[, -1L] - y, 2L, list), unlist))),
          w)),
      K <- t(backsolve(t(P.y), forwardsolve(P.y, t(P.x.y)))),
      x <- x + drop(K %*% (z - y)),
      P <- Reduce(
        ukf.chdd,
        lapply(apply(K %*% P.y, 2L, list), unlist), P),
      list(x, P)))
  
  return(xP)
}

.aux.envir.deploy <- function(env) {
  
  assign(
    "origin",
    proc.time()["elapsed"],
    envir = env)
  
  delayedAssign(
    "msg",
    c("complete\n",
      "A-GPS State Estimator\n",
      "Method to be used: %s\n",
      "Iteration %0*d/%0*d ",
      "Elapsed time: %d:%02d:%02d\n",
      "Grabbed size: %d bytes\n"),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "msg.box",
    switch(
      1L + output,
      as.function(
        alist(
          ... = ,
          invisible(NULL))),
      as.function(
        alist(
          num = , complete = , ... = ,
          `{`(
            cat(
              switch(
                1L + is.null(complete),
                paste(
                  paste(
                    rep(".", 45L - nchar(sprintf(msg[num], ...))),
                    collapse = ""),
                  msg[complete]),
                sprintf(msg[num], ...)),
              sep = ""),
            flush.console())))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "adjust",
    as.function(
      alist(
        v = ,
        k = integer(),
        switch(
          1L + is.list(k),
          replace(v, k, v[k] / sqrt(sum(v[k] ^ 2L))),
          Map(sys.function(0L), v, k)))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "schud",
    as.function(
      alist(
        S = , x = ,
        sapply(
          seq_along(x),
          as.function(
            alist(
              i = ,
              `{`(
                k <- i + seq_along(x[-1L]),
                w <- sqrt(S[i, i] ^ 2L + x[1L] ^ 2L),
                cs <- c(w, x[1L]) / S[i, i],
                out <- c(numeric(i - 1L), w,
                         (S[k, i] + cs[2L] * x[-1L]) / cs[1L]),
                x <<- cs[1L] * x[-1L] - cs[2L] * out[k],
                out)))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "schdd",
    as.function(
      alist(
        S = , x = ,
        `{`(
          y <- forwardsolve(S, x),
          norm.y <- sqrt(sum(y ^ 2L)),
          stopifnot(norm.y <= 1),
          alpha <- sqrt(1 - norm.y ^ 2L),
          cs <- t(
            apply(
              sapply(
                rev(y),
                as.function(
                  alist(
                    v = ,
                    `{`(
                      scale <- alpha + abs(v),
                      ab <- c(alpha, v) / scale,
                      n <- sqrt(sum(ab ^ 2L)),
                      alpha <<- scale * n,
                      ab / n)))),
              1L, rev)),
          t(
            sapply(
              seq_along(y),
              as.function(
                alist(
                  j = ,
                  w = 0,
                  c(
                    rev(
                      sapply(
                        j : 1L,
                        as.function(
                          alist(
                            i = ,
                            `{`(
                              out <- cs[1L, i] * S[j, i] - cs[2L, i] * w,
                              w <<- cs[1L, i] * w + cs[2L, i] * S[j, i],
                              out))))),
                    numeric(length(y) - j))))))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "rk4",
    as.function(
      alist(
        x = , tau = , fun = , ... = ,
        `+`(
          `*`(
            Reduce(
              `+`,
              Map(
                `*`,
                Reduce(
                  as.function(
                    alist(
                      p = , q = ,
                      fun(x + q * p, ...))),
                  tau * c(0, .5, .5, 1), 0, FALSE, TRUE),
                c(0, 1, 2, 2, 1))),
            tau / 6),
          x))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
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
    eval.env = env, assign.env = env)
  
  delayedAssign(
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
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "x.dot.fun",
    as.function(
      alist(
        x = , u = , w = ,
        c(x[4L : 6L],
          `+`(drop(t(S.fun(x[7L : 10L])) %*% (u[4L : 6L] + w[4L : 6L])),
              c(0, 0, params$g)),
          drop(`%*%`(.5 * Omega.fun(x[7L : 10L]),
                     u[1L : 3L] + w[1L : 3L]))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "x.fun",
    as.function(
      alist(
        tau = ,
        as.function(
          alist(
            x = , u = , w = ,
            adjust(rk4(x, tau, x.dot.fun, u, w), 7L : 10L))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "magn.x.fun",
    as.function(
      alist(
        tau = ,
        as.function(
          alist(
            x = , u = , w = ,
            x + w)))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "z.fun",
    as.function(
      alist(
        tau = , avail = ,
        as.function(
          alist(
            x = , v = , u = ,
            `+`(
              `[`(
                c(x[1L : 6L],
                  drop(S.fun(x[7L : 10L]) %*% params$E),
                  -x[3L]),
                avail),
              v))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "magn.z.fun.auto",
    as.function(
      alist(
        x = , u = ,
        sum(
          `^`(
            drop(
              `%*%`(
                matrix(
                  c(x[1L : 3L], x[2L], x[4L : 5L], x[3L], x[5L : 6L]),
                  3L , 3L),
                u[1L : 3L] - x[7L : 9L])),
            2L)))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "magn.z.fun.aided",
    as.function(
      alist(
        x = , u = ,
        `+`(
          drop(
            Reduce(
              `%*%`,
              list(
                matrix(
                  c(x[1L : 3L], x[2L], x[4L : 5L], x[3L], x[5L : 6L]),
                  3L , 3L),
                S.fun(u[4L : 7L]),
                u[1L : 3L]))),
          x[7L : 9L]))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "magn.z.fun",
    as.function(
      alist(
        tau = , avail = ,
        as.function(
          alist(
            x = , v = , u = ,
            `+`(
              `[`(
                c(magn.z.fun.auto(x, v, u),
                  magn.z.fun.aided(x, v, u)),
                avail),
              v))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.0.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (+x[07L] * u[4L] - x[10L] * u[5L] + x[09L] * u[6L]))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.1.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (+x[08L] * u[4L] + x[09L] * u[5L] + x[10L] * u[6L]))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.2.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (-x[09L] * u[4L] + x[08L] * u[5L] + x[07L] * u[6L]))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.3.fun",
    as.function(
      alist(
        x = , u = ,
        2 * (-x[10L] * u[4L] - x[07L] * u[5L] + x[08L] * u[6L]))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.V.fun",
    as.function(
      alist(
        x = , u = ,
        matrix(
          c(+J.0.fun(x, u), +J.1.fun(x, u),
            +J.2.fun(x, u), +J.3.fun(x, u),
              -J.3.fun(x, u), -J.2.fun(x, u),
              +J.1.fun(x, u), +J.0.fun(x, u),
                +J.2.fun(x, u), -J.3.fun(x, u),
                -J.0.fun(x, u), +J.1.fun(x, u)),
          3L, 4L, TRUE))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.q.fun",
    as.function(
      alist(
        x = , u = ,
        .5 * matrix(
          c(0.0E+0, -u[1L], -u[2L], -u[3L],
            +u[1L], 0.0E+0, +u[3L], -u[2L],
            +u[2L], -u[3L], 0.0E+0, +u[1L],
            +u[3L], +u[2L], -u[1L], 0.0E+0),
          4L, 4L, TRUE))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "J.fun",
    as.function(
      alist(
        tau = ,
        as.function(
          alist(
            x = , u = ,
            `*`(rbind(cbind(matrix(0, 3L, 3L),
                            diag(3L),
                            matrix(0, 3L, 4L)),
                      cbind(matrix(0, 7L, 6L),
                            rbind(J.V.fun(x, u),
                                  J.q.fun(x, u)))),
                tau))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "G.fun",
    as.function(
      alist(
        tau = ,
        as.function(
          alist(
            x = , u = ,
            `*`(cbind(rbind(matrix(0, 6L, 3L),
                            .5 * Omega.fun(x[7L : 10L])),
                      rbind(matrix(0, 3L, 3L),
                            t(S.fun(x[7L : 10L])),
                            matrix(0, 4L, 3L))),
                tau))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.0.fun",
    as.function(
      alist(
        x = ,
        `*`(2,
            `+`(`+`(+x[07L] * params$E[1L],
                    +x[10L] * params$E[2L]),
                -x[09L] * params$E[3L])))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.1.fun",
    as.function(
      alist(
        x = ,
        `*`(2,
            `+`(`+`(+x[08L] * params$E[1L],
                    +x[09L] * params$E[2L]),
                +x[10L] * params$E[3L])))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.2.fun",
    as.function(
      alist(
        x = ,
        `*`(2,
            `+`(`+`(-x[09L] * params$E[1L],
                    +x[08L] * params$E[2L]),
                -x[07L] * params$E[3L])))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.3.fun",
    as.function(
      alist(
        x = ,
        `*`(2,
            `+`(`+`(-x[10L] * params$E[1L],
                    +x[07L] * params$E[2L]),
                +x[08L] * params$E[3L])))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.B.fun",
    as.function(
      alist(
        x = ,
        matrix(
          c(+H.0.fun(x), +H.1.fun(x),
            +H.2.fun(x), +H.3.fun(x),
              +H.3.fun(x), -H.2.fun(x),
              +H.1.fun(x), -H.0.fun(x),
                -H.2.fun(x), -H.3.fun(x),
                +H.0.fun(x), +H.1.fun(x)),
          3L, 4L, TRUE))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.fun",
    as.function(
      alist(
        tau = , avail = ,
        as.function(
          alist(
            x = ,
            `[`(
              cbind(rbind(diag(6L),
                          matrix(0, 3L, 6L),
                          c(rep(0, 2L), -1, rep(0, 3L))),
                    rbind(matrix(0, 6L, 4L),
                          H.B.fun(x),
                          rep(0, 4L))),
              avail, ))))),
    eval.env = env, assign.env = env)
}

agps.estimate <- function(data,
                          sd     = rep(.1, 26L),
                          x.0    = as.numeric(data[1L, 2L : 11L]),
                          params = list(g = 9.80665,
                                        E = c(+1.94511E+2,
                                              -2.83000E+0,
                                              +4.46352E+2)),
                          method = ekf,
                          output = FALSE,
                          ...) {
  
  .aux.envir.deploy(environment())
  
  N <- length(data$time[-1L])
  fmtlen <- 1L + as.integer(trunc(log10(N)))
  msg.box(2L, NULL)
  msg.box(3L, NULL, deparse(substitute(method)))
  
  P <- sd[01L : 10L]
  Q <- sd[11L : 16L]
  R <- sd[17L : 26L]
  
  #B.sq <- sum(E ^ 2L)
  #B.para <- list() # INSERT HERE NEW ADDITIONALLY INCLUDED PORTION OF ERROR VARIANCES
  
  x <- do.call(
    rbind,
    Reduce(
      as.function(
        alist(
          v = , i = ,
          `{`(
            msg.box(4L, NULL, fmtlen, i, fmtlen, N),
            tau <- diff(data$time[i + 0L : 1L]),
            u <- as.numeric(data[i, 12L : 17L]),
            z <- as.numeric(data[i, c(2L : 7L, 18L : 21L)]),
            avail <- which(!is.na(z)),
            res <- adjust(
              method(
                v, P, Q, R[avail], u, z[avail], x.fun(tau), z.fun(tau, avail),
                ekf.J.fun = J.fun(tau),
                ekf.G.fun = G.fun(tau),
                ekf.H.fun = H.fun(tau, avail),
                ukf.chud  = schud,
                ukf.chdd  = schdd,
                ...),
              list(7L : 10L, integer())),
            P <<- res[[2L]],
            msg.box(4L, 1L, fmtlen, i, fmtlen, N),
            res[[1L]]))),
      seq_along(data$time[-1L]), x.0, FALSE, TRUE))
  
  forceAndCall(
    1L,
    as.function(
      alist(
        d = ,
        msg.box(5L, NULL, d %/% 3600L, (d %/% 60L) %% 60L, d %% 60L))),
    as.integer(round(proc.time()["elapsed"] - origin)))
  
  msg.box(6L, NULL, object.size(x))
  
  return(data.frame(time = data$time,
                    P    = x[, 01L : 03L],
                    V    = x[, 04L : 06L],
                    q    = x[, 07L : 10L]))
}

agps.simulate <- function(magnitude = list(rep(1, 3L),
                                           c(1, 1, .1)),
                          frequency = list(c(2 * pi, 2 * pi, 2 * pi),
                                           c(pi / 5, pi / 5, pi)),
                          phase     = list(rep(0, 3L),
                                           c(0, pi / 2, 0)),
                          sd        = rep(0, 16L),
                          q.0       = c(1, 0, 0, 0),
                          duration  = 10,
                          n         = 1L + 1E+3L,
                          params    = list(g = 9.80665,
                                           E = c(+1.94511E+2,
                                                 -2.83000E+0,
                                                 +4.46352E+2))) {
  
  .aux.envir.deploy(environment())
  
  omega.fun <- as.function(
    alist(
      t = ,
      magnitude[[1L]] * cos(drop(frequency[[1L]] %o% t) + phase[[1L]])))
  
  q.dot.fun <- as.function(
    alist(
      q = , u = ,
      drop(.5 * Omega.fun(q) %*% u)))
  
  time <- seq(0, duration, length.out = n)
  omega <- omega.fun(time)
  q <- cbind(q.0,
             sapply(
               seq_len(n - 1L),
               as.function(
                 alist(
                   k = ,
                   adjust(
                     rk4(q.0,
                         diff(time[k + 0L : 1L]),
                         q.dot.fun,
                         omega[, k]),
                     1L : 4L) ->> q.0))),
             deparse.level = 0L)
  S <- lapply(seq_len(n), as.function(alist(k = , S.fun(q[, k]))))
  P <- `*`(magnitude[[2L]],
           cos(drop(frequency[[2L]] %o% time) + phase[[2L]]))
  V <- `*`(-magnitude[[2L]] * frequency[[2L]],
           sin(drop(frequency[[2L]] %o% time) + phase[[2L]]))
  a <- mapply(`%*%`, S,
              unlist(
                apply(
                  `-`(`*`(-magnitude[[2L]] * frequency[[2L]] ^ 2L,
                          cos(drop(frequency[[2L]] %o% time) + phase[[2L]])),
                      c(0, 0, params$g)),
                  2L, list),
                FALSE))
  B <- sapply(S, `%*%`, params$E)
  noise <- replicate(n, rnorm(length(sd), 0, sd))
  
  return(data.frame(time  = time,
                    P     = t(P + noise[1L : 3L, ]),
                    V     = t(V + noise[4L : 6L, ]),
                    q     = t(q),
                    omega = t(omega + noise[7L : 9L, ]),
                    a     = t(a + noise[10L : 12L, ]),
                    B     = t(B + noise[13L : 15L, ]),
                    A     = -P[3L, ] + noise[16L, ]))
}

## Typical usage
## Load library for 3D plotting:
# library(rgl)
## Generate sample data, perform estimation and display result as 3D plot of the position (plus calculate elapsed time):
# system.time({x <- agps.estimate(dat <- agps.simulate(sd = rep(.1, 16L))); plot3d(x$P.1, x$P.2, x$P.3, col = "navy", size = 3L, zlim = c(-1, 1), xlab = "", ylab = "", zlab = "")})
## Add initial data to the plot:
# plot3d(dat$P.1, dat$P.2, dat$P.3, col = "red", size = 3L, add = TRUE)
## Plot component-wise error for quaternion (actual minus estimated):
# matplot(dat$time, dat[, 8L : 11L] - x[, 8L : 11L], type = "l", xlab = expression(t), ylab = expression(q))
## Draw legend for previous plot (optional):
## legend("topright", c(expression(q[1]), expression(q[2]), expression(q[3]), expression(q[4])), lty = 1L : 4L, col = 1L : 4)
## Plot values of the angle between actual quaternion and estimated one:
# plot(dat$time, acos(mapply(`%*%`, unlist(apply(x[, 8L : 11L], 1L, list), FALSE), unlist(apply(dat[, 8L : 11L], 1L, list), FALSE)) / apply(x[, 8L : 11L], 1L, as.function(alist(v = , norm(as.matrix(v), "F")))) / apply(dat[, 8L : 11L], 1L, as.function(alist(v = , norm(as.matrix(v), "F"))))), type = "l", xlab = "t", ylab = expression(varphi))

# Alternative formal signature for agps.simulate() (provides more explicit than
#   default options for omega evolution, powered by N.N.):
# agps.simulate <- function(magnitude = list(c(.1, 1, 1 / 25),
#                                            c(1, 1, .1)),
#                           frequency = list(-c(2 * pi, 2 * pi / 10, 0),
#                                            c(pi / 5, pi / 5, pi)),
#                           phase     = list(c(pi / 2, pi / 2, 0),
#                                            c(0, pi / 2, 0)),
#                           sd        = rep(0, 16L),
#                           q.0       = c(1, 0, 0, 0),
#                           duration  = 10,
#                           n         = 1L + 1E+3L,
#                           params    = list(g = 9.80665,
#                                            E = c(1.93254E-2,
#                                                  2.45820E-3,
#                                                  4.64090E-2)))
# Substitute default interface with representation specified above.

# Interesting Note. To make a function from name use, for example, the following:
#   eval(parse(text = "sin"))

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
# 
# 
# .schud.rec <- function(S, x, a = 1) {
#   
#   S <- switch(
#     1L + (length(x) != 0L), S,
#     `{`(
#       k <- seq(nrow(S) - length(x) + 1L, along.with = x),
#       w <- sqrt(S[k[1L], k[1L]] ^ 2L + a * x[1L] ^ 2L),
#       cs <- c(w, x[1L]) / S[k[1L], k[1L]],
#       S[k[1L], k[1L]] <- w,
#       S[k[-1L], k[1L]] <- (S[k[-1L], k[1L]] + a * cs[2L] * x[-1L]) / cs[1L],
#       x <- cs[1L] * x[-1L] - cs[2L] * S[k[-1L], k[1L]],
#       sys.function(0L)(S, x, a)))
#   
#   return(S)
# }
# 
# .schdd.rec <- function(S, x) {
#   
#   y <- forwardsolve(S, x)
#   norm.y <- sqrt(sum(y ^ 2L))
#   stopifnot(norm.y <= 1)
#   alpha <- sqrt(1 - norm.y ^ 2L)
#   
#   cs <- as.function(
#     alist(
#       z = rev(y),
#       aleph = alpha,
#       out = numeric(),
#       switch(
#         1L + (length(z) != 0L), out,
#         `{`(
#           scale <- aleph + abs(z[1L]),
#           ab <- c(aleph, z[1L]) / scale,
#           n <- sqrt(sum(ab ^ 2L)),
#           sys.function(0L)(z[-1L], scale * n, cbind(ab / n, out))))))()
#   
#   S <- as.function(
#     alist(
#       Q = S,
#       j = seq_along(y),
#       switch(
#         1L + (length(j) != 0L), Q,
#         sys.function(0L)(
#           as.function(
#             alist(
#               P = Q,
#               i = j[1L] : 1L,
#               w = 0,
#               switch(
#                 1L + (length(i) != 0L), P,
#                 `{`(
#                   w.new <- cs[1L, i[1L]] * w + cs[2L, i[1L]] * P[j[1L], i[1L]],
#                   P[j[1L], i[1L]] <- cs[1L, i[1L]] * P[j[1L], i[1L]] - cs[2L, i[1L]] * w,
#                   sys.function(0L)(P, i[-1L], w.new)))))(),
#           j[-1L]))))()
#   
#   return(S)
# }
#
#
# LEGACY for / oldfor COMBINATION:
#
#   assign("oldfor", `for`,
#          envir = env)
#   
#   assign(
#     "for",
#     switch(
#       1L + env$output,
#       env$oldfor,
#       as.function(
#         alist(
#           loop.var = , loop.range = , loop.expr = ,
#           `{`(
#             N <- max(loop.range),
#             fmtlen <- 1L + as.integer(trunc(log10(N))),
#             env$msg.box(2L, NULL),
#             env$msg.box(3L, NULL, "env$method.name"),
#             eval(
#               call(
#                 "oldfor",
#                 substitute(loop.var),
#                 loop.range,
#                 bquote(
#                   `{`(
#                     msg.box(4L, NULL, .(fmtlen), .(substitute(loop.var)),
#                             .(fmtlen), .(N)),
#                     .(substitute(loop.expr)),
#                     msg.box(4L, 1L, .(fmtlen), .(substitute(loop.var)),
#                             .(fmtlen), .(N))))),
#               envir = env),
#             as.function(
#               alist(
#                 d = as.integer(round(proc.time()["elapsed"] - env$origin)),
#                 env$msg.box(5L, NULL,
#                             d %/% 3600L, (d %/% 60L) %% 60L, d %% 60L)))(),
#             env$msg.box(6L, NULL, object.size(env$x)))))),
#     envir = env)
# 
# EOF
# THANK YOU FOR YOUR ATTENTION!
