#' locality preserving projections
#'
#' @param X the data matrix
#' @param W the affinity matrix
#' @param ndim number of dimensions
#' @importFrom RSpectra eigs
#' @export
embed_lpp <- function(X, W, ndim=2) {

  D <- Diagonal(x=rowSums(W))
  L <- D - W

  Zl <- t(X) %*% L %*% X
  Zr <- t(X) %*% D %*% X

  decomp <- RSpectra::eigs(solve(Zr,Zl), k=ndim, which="SM")

}
