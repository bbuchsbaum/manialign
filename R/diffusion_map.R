


diffusion_map <- function(A, ncomp=2) {

  n <- nrow(A)

  D <- rowSums(A)
  Ds <- 1/Diagonal(x=sqrt(D))

  ## correct?
  L <- Diagonal(x=rep(1,n)) - Ds %*% A %*% Ds

  # make A matrix sparse
  #ind = which(A>delta, arr.ind=TRUE)
  #Asp = sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims = c(n,n))

  decomp <- RSpectra::svds(L, k=ncomp)
  cat('Performing eigendecomposition\n') # eigendecomposition


  # eigendecomposition using ARPACK
  decomp = arpack(f,extra=Asp,sym=TRUE,
                options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
  psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
  phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
  eigenvals = decomp$values #eigenvalues


cat('Computing Diffusion Coordinates\n')
if(t<=0){# use multi-scale geometry
  lambda=eigenvals[-1]/(1-eigenvals[-1])
  lambda=rep(1,n)%*%t(lambda)
  if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
    lam = lambda[1,]/lambda[1,1]
    neigen = min(which(lam<.05)) # default number of eigenvalues
    neigen = min(neigen,maxdim)
    eigenvals = eigenvals[1:(neigen+1)]
    cat('Used default value:',neigen,'dimensions\n')
  }
  X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
}
else{# use fixed scale t
  lambda=eigenvals[-1]^t
  lambda=rep(1,n)%*%t(lambda)

  if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
    lam = lambda[1,]/lambda[1,1]
    neigen = min(which(lam<.05)) # default number of eigenvalues
    neigen = min(neigen,maxdim)
    eigenvals = eigenvals[1:(neigen+1)]
    cat('Used default value:',neigen,'dimensions\n')
  }

  X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
}
cat('Elapsed time:',signif(proc.time()[3]-start,digits=4),'seconds\n')

y = list(X=X,phi0=phi[,1],eigenvals=eigenvals[-1],eigenmult=lambda[1,1:neigen],
         psi=psi,phi=phi,neigen=neigen,epsilon=eps.val)
class(y) = "dmap"
return(y)
