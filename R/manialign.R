
## https://math.stackexchange.com/questions/1106343/generalized-eigenvalue-problem-for-symmetric-low-rank-matrix


geig <- function(A, B, ncomp=min(3,dim(A)), which="LA") {

  Xprime <- Diagonal(x=1/sqrt(diag(B)))
  #C  <- t(Xprime) %*% A %*% Xprime
  CC <- crossprod(Xprime, A) %*% Xprime
  eres <- eigs_sym(CC,symmetric=TRUE, k=ncomp, which=which)
  XX <- Xprime %*% eres$vectors
  list(vectors=XX, values=eres$values, ncomp=ncomp)
}




#' pairwise_id_matrix
#'
#'
#' @export
pairwise_label_matrix <- function(Ls, offsets, type="s") {

  out <- do.call(rbind, lapply(1:length(Ls), function(i) {
    do.call(rbind, lapply(1:length(Ls), function(j) {
      m <- neiborweights:::label_matrix(Ls[[i]], Ls[[j]], type=type, return_matrix=FALSE)
      m[,1] <- m[,1] + offsets[i]
      m[,2] <- m[,2] + offsets[j]
      m
    }))
  }))

  len <- sum(sapply(Ls, length))

  sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(len, len))
}


#' pairwise_id_matrix
#'
#'
#' @export
pairwise_id_matrix <- function(Ids, offsets, lens) {
  out <- do.call(rbind, lapply(1:length(Ids), function(i) {
    do.call(rbind, lapply(1:length(Ids), function(j) {
      if (i == j) {
        NULL
      } else {
        m <- neighborweights:::label_matrix(Ids[[i]], Ids[[j]],
                                            type="s", return_matrix=FALSE,
                                            dim1=lens[i], dim2=lens[j])
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
predict.mani_align_correspondence <- function(x, newdata, table_index, ncomp=x$ncomp) {
  ncomp <- min(x$ncomp, ncomp)
  ind <- x$block_indices[[table_index]]
  proj <- x$vectors[ind, 1:ncomp]
  t(proj) %*% t(newdata)
}


gen_correspondence_laplacian <- function(Xs, id_set, k=10, sigma=.73) {

  ninstances <- unlist(lapply(Xs, ncol))
  nsets <- length(Xs)
  offsets <- cumsum(c(0, ninstances[1:(nsets-1)]))
  message("setting up id matrix")
  Wc <- pairwise_id_matrix(id_set, offsets, lens=ninstances)
  message("computing knn")
  Ws <- Matrix::bdiag(lapply(Xs, function(x) {
    message("sim mat")
    neighborweights::similarity_matrix(t(x),neighbor_mode="knn",k=knn, sigma=sigma)
  }))

  Ds <- Diagonal(x=rowSums(Ws))
  #Di <- Diagonal(x=rep(1, nrow(Ws)))
  Omega <- Diagonal(x=rowSums(Wc))

  list(Omega=Omega, Ds=Ds, Ws=Ws, Wc=Wc)

}

#' @export
#' @param Xs a list of matrices to be aligned, where rows are instances and columns features.
#' @param id_set a list of ids for the columns (instances) of each matrix used to derive correspondences.
#' @param ncomp the number of components to estimate.
#' @param knn the number nearest neighbors used to define the similarity matrices.
#' @param sigma the bandwidth of the heat kernel
#' @param u1 the weight placed on topology preservation.
#' @param u2 the weight placed on alignment.
#' @importFrom neighborweights similarity_matrix
mani_align_instances <- function(Xs, id_set, ncomp=2, knn=10, sigma=.73, u1=.5, u2=.5) {
  C <- gen_correspondence_laplacian(Xs, id_set, k=knn, sigma)
  Ls <- C$Ds - C$Ws
  L <- u2*Ls + u1*C$Omega -u1*C$Wc

  geig(L, C$Ds, ncomp=ncomp)
  #Matrix::nearPD(L)
  decomp <- RSpectra::eigs(solve(C$Ds, L), ncomp)
  ret <- list(vectors=decomp$vectors,
              values=decomp$values,
              block_indices=block_indices,
              knn=k,
              sigma=sigma,
              u1=u1,
              u2=u2,
              ncomp=ncomp)

  class(ret) <- c("mani_align_instances", "list")
  ret
}


#' @export
#' @importFrom neighborweights similarity_matrix
mani_align_features <- function(Xs, id_set, ncomp=2, knn=10, sigma=.73, u1=.5, u2=.5) {
  C <- gen_correspondence_laplacian(Xs, id_set, sigma)


  Ls <- C$Ds - C$Ws
  L <- u2*Ls + u1*C$Omega -u1*C$Wc

  Z <- Matrix::bdiag(Xs)

  Zl <- Z %*% L %*% t(Z)
  Zr <- Z %*% C$Ds %*% t(Z)

  decomp <- eigen(solve(Matrix::nearPD(Zr)$mat, Matrix::nearPD(Zl)$mat ))

  ret <- list(vectors=decomp$vectors,
              values=decomp$values,
              block_indices=block_indices,
              knn=k,
              sigma=sigma,
              u1=u1,
              u2=u2,
              ncomp=ncomp)

  class(ret) <- c("mani_align_features", "list")
  ret

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

  W <- bdiag(lapply(Xs, function(x) neighborweights::construct_weight_matrix(x,
                                                                             neighbor_mode="knn",
                                                                             k=knn, sigma=sigma)))

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
