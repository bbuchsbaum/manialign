#' locality preserving projections
#'
#' @param X the data matrix
#' @param W the affinity matrix
#' @param ncomp number of latent components
#' @importFrom RSpectra eigs
#' @export
lpp <- function(X, W, ncomp=2, reduced_rank=ncol(X), center=TRUE, scale=FALSE) {
  assert_that(nrow(X) == nrow(W))

  X <- neuroca::pre_processor(X)

  d <- rowSums(W)
  D <- Diagonal(x=d)
  D_mhalf <- Diagonal(x=d^(-1/2))

  W <- D_mhalf %*% W %*% D_mhalf
  #W <- neighborweights:::psparse(W, max)

  Xp <- D^.5 %*% X

  #L <- D - W
  #Zl <- t(X) %*% L %*% X
  #Zr <- t(X) %*% D %*% X
  #decomp <- RSpectra::eigs(solve(Zr,Zl), k=ndim, which="SM")

  res <- lge(Xp, W, D=NULL, ncomp=ncomp, reduced_rank=reduced_rank)

}
