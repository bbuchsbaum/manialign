psparse <- function(..., fun=c("max", "min"), na.rm=FALSE) {
  fun <- match.arg(fun)
  # check that all matrices have conforming sizes
  num.rows <- unique(sapply(list(...), nrow))
  num.cols <- unique(sapply(list(...), ncol))
  stopifnot(length(num.rows) == 1)
  stopifnot(length(num.cols) == 1)

  cat.summary <- do.call(rbind, lapply(list(...), summary))


  out.summary <- if (fun == "min") {
    aggregate(x ~ i + j, data = cat.summary, FUN=function(x) {
      if (length(x) == 1) 0 else x[1]
    })
  } else {
    aggregate(x ~ i + j, data = cat.summary, FUN=max, na.rm)
  }

  sparseMatrix(i = out.summary$i,
               j = out.summary$j,
               x = out.summary$x,
               dims = c(num.rows, num.cols))
}

pmin.sparse <- function(..., na.rm = FALSE) { psparse(..., fun="min") }

pmax.sparse <- function(..., na.rm = FALSE) { psparse(..., fun="max" ) }


#' @export
sim_knn_from_adj <- function(A, k=5, type=c("normal", "mutual"), ncores=1) {
  assert_that(k > 0 && k <= nrow(X))

  type <- match.arg(type)

  jind <- 1:nrow(X)
  A2 <- do.call(rbind, mclapply(1:nrow(A), function(i) {
    ord <- order(A[i,], decreasing=TRUE)
    cbind(i=i, j=jind[ord[1:k]],x=A[i, ord[1:k]])
  }, mc.cores=ncores))

  m <- sparseMatrix(i=A2[,1], j=A2[,2], x=A2[,3], dims=c(nrow(A), nrow(A)))
  m <- if (type == "normal") {
    pmax.sparse(m, t(m))
  } else {
    pmin.sparse(m, t(m))
  }
}

#' @importFrom assertthat assert_that
#' @importFrom FNN get.knnx
#' @importFrom parallel mclapply
#' @export
sim_knn <- function(X, k=5, sigma=1, type=c("normal", "mutual"), ncores=1) {
  assert_that(k > 0 && k <= nrow(X))

  type <- match.arg(type)
  nn <- get.knnx(X, X, k=k+1)

  jind <- 1:nrow(X)
  A <- do.call(rbind, mclapply(1:nrow(X), function(i) {
    cbind(i=i, j=jind[nn$nn.index[i,2:k]],x=exp(-nn$nn.dist[i,2:k]/sigma^2))
  }, mc.cores=ncores))

  m <- sparseMatrix(i=A[,1], j=A[,2], x=A[,3], dims=c(nrow(X), nrow(X)))
  if (type == "normal") {
    pmax.sparse(m, t(m))
  } else {
    pmin.sparse(m, t(m))
  }
}
