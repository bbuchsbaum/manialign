





#' @importFrom Matrix sparseMatrix
label_matrix <- function(a, b, type=c("s", "d"), return_matrix=TRUE, simfun=NULL , dim1=length(a), dim2=length(b)) {
  type <- match.arg(type)

  a.idx <- which(!is.na(a))
  b.idx <- which(!is.na(b))

  sfun <- if (!is.null(simfun)) {
    simfun
  } else {
    function(x, y) {
      if (is.na(x) || is.na(y)) {
        0
      } else if (x == y) {
        1
      } else {
        0
      }
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
    sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(dim1, dim2))
  } else {
    out
  }
}


pairwise_label_matrix <- function(Ls, offsets, type="s") {

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


pairwise_id_matrix <- function(Ids, offsets, lens) {
  out <- do.call(rbind, lapply(1:length(Ids), function(i) {
    do.call(rbind, lapply(1:length(Ids), function(j) {
      if (i == j) {
        NULL
      } else {
        m <- label_matrix(Ids[[i]], Ids[[j]], type="s", return_matrix=FALSE, dim1=lens[i], dim2=lens[j])
        m[,1] <- m[,1] + offsets[i]
        m[,2] <- m[,2] + offsets[j]
        m
      }
    }))
  }))

  len <- sum(lens)
  sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(len, len))

}


#' @export
#' @importFrom neighborweights construct_weight_matrix
mani_align_correspondence <- function(Xs, id_set, ncomp=2, knn=10, sigma=.73, u1=.5, u2=.5) {
  ninstances <- unlist(lapply(Xs, nrow))
  nsets <- length(Xs)

  offsets <- cumsum(c(0, ninstances[1:(nsets-1)]))
  block_indices <- sapply(1:length(ninstances), function(i) {
    seq(offsets[i] + 1, offsets[i] + ninstances[i])
  })

  Wc <- pairwise_id_matrix(id_set, offsets, lens=ninstances)
  Ws <- Matrix::bdiag(lapply(Xs, function(x) neighborweights::construct_weight_matrix(x, neighbor_mode="knn", k=knn, sigma=sigma)))



  Ds <- Diagonal(x=rowSums(Ws))
  Omega <- Diagonal(x=rowSums(Wc))


  Ls <- Ds - Ws
  L <- u2*Ls + u1*Omega -u1*Wc

  Z <- Matrix::bdiag(Xs)

  eigen(solve(B,A) )$vectors

  Zl <- Z %*% L %*% t(Z)
  Zr <- Z %*% Ds %*% t(Z)

  decomp <- eigen(solve(Zr,Zl) )
  ret <- list(vectors=decomp$vectors,
              values=decomp$values,
              block_indices=block_indices,
              knn=k,
              sigma=sigma,
              u1=u1,
              u2=u2)

  class(ret) <- c("mani_align_correspondence", "list")


}


#' @export
mani_align_labeled <- function(Xs, label_set, ncomp=2, knn=5, sigma=.73, u1=1, u2=1) {
  ninstances <- unlist(lapply(Xs, nrow))
  nsets <- length(Xs)

  offsets <- cumsum(c(0, ninstances[1:(nsets-1)]))

  block_indices <- sapply(1:length(ninstances), function(i) {
    seq(offsets[i] + 1, offsets[i] + ninstances[i])
  })

  Ws <- pairwise_label_matrix(label_set, offsets, type="s")
  Wd <- pairwise_label_matrix(label_set, offsets, type="d")

  neighborweights::construct_weight_matrix(x, neighbor_mode="knn", k=knn, sigma=sigma)
  W <- bdiag(lapply(Xs, function(x) neighborweights::construct_weight_matrix(x, neighbor_mode="knn", k=knn, sigma=sigma)))

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
