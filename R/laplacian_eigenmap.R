
#' @importFrom Matrix rowSums
eigenmap <- function(A, normalized=TRUE, ncomp=2) {
  stopifnot(isSymmetric(A))
  stopifnot(ncomp > 0 & ncomp < nrow(A))

  D <- rowSums(A)
  L <- Diagonal(x=D) - A

  decomp <- if (normalized) {
    Dtilde <- Diagonal(x= D^(-1/2))
    NL <- Ds%*%L%*%Ds
    RSpectra::eigs(NL, which="SM", k=ncomp)
  } else {
    geig(L, D, which="SM", ncomp=ncomp)
  }

  nc <- length(decomp$values)
  ret <- list(vectors=decomp$vectors[,1:(nc-1)], values=decomp$values[1:nc-1])
  class(ret) <- c("eigenmap", "list")
  ret
}




