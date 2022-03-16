mag <- function()
{
  raw <- read.table("NDU_RAW_data.txt")
  B <- unname(as.matrix(raw[, 13L : 15L]))
  X <- cbind(B, 1)
  Y <- apply(B ^ 2L, 1L, sum)
  R <- drop(solve(t(X) %*% X) %*% t(X) %*% Y)
  return(list(V = R[1L : 3L] / 2, B = sqrt(R[4L] + sum((R[1L : 3L] / 2) ^ 2L))))
}

mag2 <- function()
{
  raw <- read.table("NDU_RAW_data.txt")
  X <- unname(as.matrix(raw[, 13L : 15L]))
  A <- function(x)
    do.call(rbind, list(x[1L : 3L], x[c(2L, 4L : 5L)], x[c(3L, 5L : 6L)]))
  dA <- function(x)
    do.call(c, lapply(1L : length(x), as.function(alist(i = , x[i] * x[i : length(x)]))))
  G <- function(x)
    apply(X, 1L, as.function(alist(z = , drop((z - x[1L : 3L]) %*% A(x[4L : 9L]) %*% (z - x[1L : 3L])) - x[10L] ^ 2L)))
  J <- function(x)
    t(
      apply(
        X,
        1L,
        as.function(
          alist(
            z = ,
            c(
              A(x[4L : 9L]) %*% (z - x[1L : 3L]) + t(A(x[4L : 9L])) %*% (z - x[1L : 3L]),
              dA(z - x[1L : 3L]),
              - 2 * x[10L])))))
  x.0 <- mag()
  x <- c(x.0$V, c(1, 0, 0, 1, 0, 1), x.0$B)
  gr <- numeric(0L)
  G.c <- G(x)
  for(i in 1L : 21L)
  {
    J.c <- J(x)
    x <- x - chol2inv(chol(t(J.c) %*% J.c)) %*% t(J.c) %*% G.c
    G.c <- G(x)
    gr <- c(gr, norm(2 * t(J.c) %*% G.c, "F"))
    cat(sprintf("%.3d", i), ":", sprintf("%+.10E", drop(G.c %*% G.c)), "/", sprintf("%+.10E", gr[length(gr)]), "\n")
    flush.console()
  }
  return(list(V = x[1L : 3L], A = A(x[4L : 9L]), B = x[10L]))
}