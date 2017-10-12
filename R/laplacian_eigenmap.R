
#' @importFrom Matrix rowSums
laplacian_eigenmap <- function(A, normalized=TRUE, ncomp=2) {
  stopifnot(isSymmetric(A))
  stopifnot(ncomp > 0 & ncomp < nrow(A))

  D <- rowSums(A)
  L <- Diagonal(x=D) - A

  decomp <- if (normalized) {
    Dtilde <- Diagonal(x= 1/(D^2))
    NL <- Ds%*%L%*%Ds
    RSpectra::eigs(NL, which="SM", k=ncomp)
  } else {
    geig(L, D, which="SM", ncomp=ncomp)
  }

  nc <- length(decomp$values)
  ret <- list(vectors=decomp$vectors[,1:(nc-1)], values=decomp$values[1:nc-1])
  class(ret) <- c("laplacian_eigenmap", "list")
  ret
}




#' @importFrom RSpectra eigs
commute_time <- function(A, ncomp=nrow(A)) {
  D <- Matrix::rowSums(A)

  Dtilde <- Diagonal(x= 1/(D^2))

  M <- Dtilde %*% A %*% Dtilde

  decomp <- RSpectra::eigs(M, k=ncomp)
  pii <- D/sum(A)
  v <- decomp$vectors[, 2:ncomp]

  cds <- sweep(v, 2, sqrt(1 - decomp$values[2:ncomp]), "/")
  cds <- sweep(cds, 2, 1/pi)

}
