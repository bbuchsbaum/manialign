


diffusion_map <- function(A, ncomp=2, t=0) {

  n <- nrow(A)

  D <- rowSums(A)
  Ds <- Diagonal(x=1/sqrt(D))
  As <- Ds %*% A %*% Ds

  ## correct?
  ## L <- Diagonal(x=rep(1,n)) - Ds %*% A %*% Ds

  # make A matrix sparse
  #ind = which(A>delta, arr.ind=TRUE)
  #Asp = sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims = c(n,n))

  ncomp <- min(ncomp, nrow(A) -1)

  decomp <- RSpectra::eigs_sym(As, k=ncomp+1)
  cat('Performing eigendecomposition\n') # eigendecomposition

  psi <- decomp$vectors / (decomp$vectors[,1]%*%matrix(1,1,ncomp+1))#right ev
  phi <- decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,ncomp+1))#left ev
  eigenvals <- decomp$values #eigenvalues
  n <- nrow(A)

  lambda <- if(t<=0){
    # use multi-scale geometry
    eigenvals[-1]/(1-eigenvals[-1])
  } else{
    # use fixed scale t
    eigenvals[-1]^t
  }

  lambda <- rep(1,n)%*%t(lambda)
  cds <- psi[,2:(ncomp+1)]*lambda[,1:ncomp] #diffusion coords. X

  ret <- list(cds=cds,phi0=phi[,1],
              eigenvals=eigenvals[-1],
              eigenmult=lambda[1,1:ncomp],
              psi=psi,phi=phi,ncomp=ncomp)
  class(ret) <- c("diffusion_map", "list")
  ret
}
