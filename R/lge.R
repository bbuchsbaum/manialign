
#' linear graph embedding
#'
#' @param X the data matrix with with rows as samples and columns as features
#' @param W the affinity matrix
#' @param D the constraint matrix
#' @param ncomp the number of eigenvpairs to retain
#' @param reduced_rank the number of PCA dimensions to retain
lge <- function(X, W, D=NULL, ncomp=3, reduced_rank=ncol(X)) {
  nsamples <- nrow(X)
  nfeatures <- ncol(X)

  if (nrow(W) != nsamples) stop('W and X mismatch!');

  if (!is.null(D) && (nrow(D) != nsamples)) {
    stop('D and data mismatch!')
  }


  if (reduced_rank < ncol(X)) {
    is_reduced <- TRUE

    pca_fit <- neuroca::pca(X, ncomp=reduced_rank, center=FALSE, scale=FALSE)

    if (!is.null(D)) {
      X <- scores(pca_fit)
      pca_proj <- loadings(pca_fit)
      DPrime <- t(X) %*% D %*% X

      ## symmetrizes DPrime
      DPrime = neighborweights:::psparse(DPrime, FUN=max)
    } else {
      X <- pca_fit$u
      pca_proj = loadings(pca_fit) %*% Matrix::Diagonal(x=1/singular_values(pca_fit))
    }
  } else {
    is_reduced <- FALSE
    pca_proj <- NULL
    DPrime <- if (!is.null(D)) {
      t(X) %*% D %*% X
    } else {
      t(X) %*% X
    }
  }

  WPrime <- t(X)  %*% W %*% X
  WPrime <- neighborweights:::psparse(WPrime,max)

  decomp <- if (is_reduced && is.null(D)) {
    RSpectra::eigs(WPrime,k=ncomp,'LA')
  } else {
    geig(WPrime, DPrime, ncomp=ncomp,'LA')
  }

  eigvector <- decomp$vectors

  if (is_reduced) {
    eigvector = pca_proj %*% eigvector
  }

  for (i in 1:ncol(eigvector)) {
    eigvector[,i] = eigvector[,i]/norm(eigvector[,i], "2")
  }

  ret <- list(vectors=eigvector, values=decomp$values, pca_proj=pca_proj)
  class(ret) <- c("lge", "list")
  ret

}

lsda <- function(X, labels, k=5, sigma=.7, beta=.1, alpha=.1, reduced_rank=ncol(X)) {
  labels <- as.factor(labels)
  nlabel <- length(levels(labels))

  Ww <- neighborweights:::label_matrix(labels, labels, type="s")
  Wb <- neighborweights:::label_matrix(labels, labels, type="d")

  if (k > 0) {
    W <- neighborweights::similarity_matrix(X, neighbor_mode = "knn", weight_mode="normalized", k=k, sigma=sigma, labels=labels)
    Ww <- Ww * W
    Wb <- Wb * W
  }

  Db <- rowSums(Wb,2)
  Wb <- -Wb

  Wb <- Wb + Matrix::Diagonal(x=Db)

  D <- rowSums(Ww)

  Wcombined = (beta/(1-beta))*Wb + Ww
  D <- Diagonal(x=D)

  res <- lge(X, Wcombined, D, reduced_rank=reduced_rank)
  class(res) <- c("lsda", "list")
  res

}

repmat <- function(a,n,m) kronecker(matrix(1,n,m),a)

