
gen_dmat <- function(nm, X) {

  S <- matrix(0, ncol(X), ncol(X))

  for (i in 1:nrow(nm)) {
    vec <- X[i,]
    idx <- which(nm[i,] > 0)
    print(length(idx))

    if (length(idx) > 0) {
      S <- S + Reduce("+", lapply(idx, function(j) {
        outer(vec-X[j,], vec-X[j,])
      }))
    }
  }

  S

}

initSw <- function(X, Y) {
  wt <- Matrix(0, nrow(X), nrow(X))
  nlevs <- length(levels(Y))
  levs <- levels(Y)

  within_class <- function(x) {
    n <- nrow(x)
    S <- n * Reduce("+", lapply(1:nrow(x), function(i) {
      Reduce("+", lapply(1:nrow(x), function(j) {
        v <- x[i,] - x[j,]
        (1/n)^2 * outer(v,v)
        #outer(v,v)
      }))
    }))
  }

  Sall <- Reduce("+", lapply(levs, function(lev) {
    idx <- which(Y == lev)
    x <- X[idx,]
    within_class(x)
  }))

}






nmmp <-function(X, Y, kw=1, kb=1, ncomp=ncol(X)-1) {
  X <- scale(X)

  #nb <- between_class_weights(X, k=3, labels=Y, weight_mode="binary")
  nw <- within_class_weights(X, k=3, labels=Y, weight_mode="binary")

  #Sw <- gen_dmat(nw, X)
  #Sb <- gen_dmat(nb, X)

  ## A = Sw
  ## B = St

  ## max St/Sw

  St <- 2*(cov(X) * nrow(X)-1)
  Sw <- initSw(X,Y)

  l1 <- matrixcalc::matrix.trace(St)/matrixcalc::matrix.trace(Sw)

  l2 <- sum(eigen(St)$values[1:ncomp])/sum(eigen(Sw)$values[(ncol(X)-ncomp+1):ncol(X)])

  lambda <- (l1 + l2)/2

  tol <- 1e-5
  while (l2 - l1 > tol) {
    ret <- eigen(St - lambda*Sw)
    print(ret$values)
    g <- sum(ret$values[1:ncomp])
    if (g > 0) {
      l1 <- lambda
    } else {
      l2 <- lambda
    }
    lambda <- (l1+l2)/2
    print(lambda)
  }

}
