


#' @importFrom Matrix sparseMatrix
label_matrix <- function(a, b, type=c("s", "d"), return_matrix=TRUE) {
  type <- match.arg(type)

  a.idx <- which(!is.na(a))
  b.idx <- which(!is.na(b))

  sfun <- function(x, y) {
    if (is.na(x) || is.na(y)) {
      0
    } else if (x == y) {
      1
    } else {
      0
    }
  }

  dfun <- function(x, y) {
    if (is.na(x) || is.na(y)) {
      0
    } else if (x == y) {
      0
    } else {
      1
    }
  }

  fun <- if (type == "s") sfun else dfun

  out <- lapply(a.idx, function(i) {
    ret <- lapply(b.idx, function(j) {
      if (fun(a[i], b[j])) {
        c(i,j,1)
      } else {
        NULL
      }
    })
    ret[sapply(ret, function(x) !is.null(x))]
  })

  out <- unlist(out, recursive=FALSE)
  out <- do.call(rbind, out)

  if (return_matrix) {
    sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(length(a), length(b)))
  } else {
    out
  }
}

#' @export
pairwise_label_matrix <- function(Ls, offsets, type) {
  out <- do.call(rbind, lapply(1:length(Ls), function(i) {
    do.call(rbind, lapply(1:length(Ls), function(j) {
      m <- label_matrix(Ls[[i]], Ls[[j]], type=type, return_matrix=FALSE)
      m[,1] <- m[,1] + offsets[i]
      m[,2] <- m[,2] + offsets[j]
      m
    }))
  }))

  len <- sum(sapply(Ls, length))

  sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(len, len))
}


#' @export
mani_align_labeled <- function(Xs, label_set, k=5, sigma=.73, ncores=1, u1=1, u2=1) {
  ninstances <- unlist(lapply(Xs, nrow))
  nsets <- length(Xs)

  offsets <- cumsum(c(0, ninstances[1:(nsets-1)]))

  block_indices <- sapply(1:length(ninstances), function(i) {
    seq(offsets[i] + 1, offsets[i] + ninstances[i])
  })

  Ws <- pairwise_label_matrix(label_set, offsets, type="s")
  Wd <- pairwise_label_matrix(label_set, offsets, type="d")
  W <- bdiag(lapply(Xs, function(x) sim_knn(x, k=k, sigma=sigma)))

  Ds <- Diagonal(x=rowSums(Ws))
  Dd <- Diagonal(x=rowSums(Wd))
  D <- Diagonal(x=rowSums(W))

  Ls <- Ds - Ws
  Ld <- Dd - Wd
  L <- D - W

  Z <- bdiag(Xs)

  Zl <- Z %*% (u1*Ls + u2*L) %*% t(Z)
  Zr <- Z %*% Ld %*% t(Z)

  decomp <- eigen(solve(Zr,Zl) )
  ret <- list(vectors=decomp$vectors,
       values=decomp$values,
       block_indices=block_indices,
       k=k,
       sigma=sigma,
       u1=u1,
       u2=u2)

  class(ret) <- c("mani_align_labels", "list")

  #B_inv <- corpcor::pseudoinverse(Zr)
  #r1 <- eigen(B_inv %*% Zl)

}
