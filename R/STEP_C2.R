# STEP C2 (generate kernel list)

# fit = fit1[[4]]; y = NULL;epsilon = 1e-4; B = "No"; verbose = FALSE; maf = 0.1;maf_beta=0.5; maf_u = 0.5

composite_kernel <- function(fit = fit, y = NULL,epsilon = 1e-4, B = "No", verbose = FALSE, maf = 0.1,
                             maf_beta=0.5, maf_u = 0.5){

  #------------------------------------------------------------------------------------------------------------
  # (0) initial

  index_test = fit$index_test # index of the NA value in Y (index for testing samples)

  # index_fixed_marker = fit$index_fixed_marker # index of markers for fixed term

  N <- fit$N # total num of samples

  M <- fit$M # num of regions in the testing data

  Px <- fit$Px # num of SNPs in the testing data

  P <- fit$DIM # num of SNPs in the each region

  w_beta <- fit$E_w_beta # weight of beta

  w_u    <- fit$E_w_u # weight of u

  # geno_fix <- fit$genotype0
  geno_fix <- fit$genotype
  # all data - fixed term

  geno_random <- lapply(1:M, function(m) attr(fit,"UD")[[m]] ) # all data - random term

  geno_background <- fit$gamma_lambda # all data - background term

  X <- fit$X  # data for training - fixed term

  Cov_List0 <- fit$Cov_List

  Z <- fit$Z # data for training - random term

  m_u <- fit$m_u # mean of u

  m_beta <- fit$m_beta # mean of beta

  m_b <- fit$m_b # mean of background

  post_sigma2e <- fit$post_sigma2e # sigma2 of error

  post_sigma2b <- fit$post_sigma2b # sigma2 of error

  var_beta <- fit$post_sigma2beta # sigma2 of error

  var_u <- fit$post_sigma2u # sigma2 of error

  yc <- fit$y

  center <- fit$center

  scale <- fit$scale

  #------------------------------------------------------------------------------------------------------------
  # (1) sorting fixed effects by w_beta

  delta_sort_PVE_beta <- rep(-Inf, Px)

  # obtain the var_beta order index

  # var_beta <- sapply(1:Px, function(j) var(as.matrix(X[,j]) %*% w_beta[j] %*% m_beta[j])  )

  # lst_beta <- sort(var_beta, index.return=TRUE, decreasing=TRUE)

  # (1) sorting random effects by w_u

  delta_sort_PVE_u <- rep(-Inf, M)

  # obtain the var_u order index

  # var_u <- sapply(1:M, function(m) var(Z[[m]] %*% (drop(w_u[[m]]) * diag(1, P[[m]])) %*% m_u[[m]])  )

  lst_u <- sort(var_u, index.return=TRUE, decreasing=TRUE)

  #------------------------------------------------------------------------------------------------------------
  # (2) fixed effects - PVE and select top markers
  # PVE - the total proportion of variance in phenotype explained by the sparse effects and random effects terms together

  # calculate and store the PVE_beta
  PVE_beta <- sapply(1:Px, function(j) var(as.matrix(X[,j]) %*% m_beta[j])/(post_sigma2e + var_beta + sum(var_u) + post_sigma2b) )

  # finial index for fixed term # select top regions based on index_delta_PVE_beta
  fix_index <- which(w_beta>=maf_beta)

  # (2) random effects - PVE and select top regions -------------------------------------------------------------
  #
  # calculate and store the PVE_u
  PVE_u <- sapply(1:M, function(m) var_u[[m]]/( post_sigma2e+var_beta + sum(var_u) + post_sigma2b) )

  #PVE_u <- sapply(1:M, function(j) var(TEST$m_u[,j])/var(y) ) # calculate and store the PVE_u

  # sort the PVE_u by the pi_u (decreasing) (add 0 to be first for further calculation)

  sort_PVE_u <- append(0, sapply(1:M, function(m) sum(PVE_u[lst_u$ix][1:m]) ))

  # index of change
  if(length(sort_PVE_u)>1){

    delta_sort_PVE_u <-  sapply(1:M, function(m) delta_sort_PVE_u[m] <- sort_PVE_u[m + 1] - sort_PVE_u[m] )

    index_delta_PVE_u <- which( abs(delta_sort_PVE_u) > epsilon )

    len <- min(length(index_delta_PVE_u), ceiling(0.8*M))

  }else if(length(sort_PVE_u)==1){

    delta_sort_PVE_u <- sort_PVE_u

    index_delta_PVE_u <- 1

    len <- 1
  }

  # finial index for random term # select top regions based on index_delta_PVE_u
  # random_index <- lapply(lst_u, `[`, lst_u$x %in% head(lst_u$x,(index_delta_PVE_u[1] - 1)) )

  random_index <- lapply(lst_u, `[`, lst_u$x %in% utils::head(lst_u$x,len) )

  #------------------------------------------------------------------------------------------------------------
  # (3) calculate the PGE - the proportion of genetic variance explained by the sparse effects terms

  # 1 get the weight for each kernel

  nr <- length(random_index$ix)

#  nr <- 1
 # Cov_List2 <- list()
  Cov_List <- list()
  wr <- list()

  if(nr!=0){
    for (i in 1:nr) {
      wr[[i]] <- var_u[random_index$ix[i]]/sum(var_u[random_index$ix])
    #  wr[[i]] <- 1
      Cov_List[[i]] <- wr[[i]] * Cov_List0[[random_index$ix[i]]]
    #  Cov_List[[i]] <- wr[[i]] * Cov_List0[[1]]
    }
    sel_ran <- which(wr>1e-2)
    composite_kernel <- Reduce('+', Cov_List[sel_ran])
   # composite_kernel <- Reduce('+', Cov_List)
  }else{
    composite_kernel <- matrix(0, 712, 712)
  }

  # 2 get composite kernel

  # 2 get additive matrix

#  mat_add <- X[,fix_index$ix]
 # mat_add <- geno_fix[,fix_index]
   mat_add <- geno_fix

   PVE_u2 <- Reduce(`+`, PVE_u) + Reduce(`+`, PVE_beta)

    return(list(mat_add = mat_add, composite_kernel = composite_kernel, wr = wr, PVE_u = PVE_u2))

}

# fit1=fit1; ng=10;

#' @importFrom utils head
STEP_C2 <- function(fit1, ng){

  ckernel <- list()
  cmat_add <- list()
  cmat_kernel <- list()
  cmat_add_dim <- list()
  wr <- list()
  PVE_u <- list()

  for (i in 1:ng) {
    ckernel[[i]] <- composite_kernel(fit = fit1[[i]], y = NULL,epsilon = 1e-4, B = "No", verbose = FALSE, maf = 0.1)
    cmat_add[[i]] <- ckernel[[i]]$mat_add
    cmat_add_dim[[i]] <- dim(cmat_add[[i]])
    cmat_kernel[[i]] <- ckernel[[i]]$composite_kernel
    wr[[i]] <- ckernel[[i]]$wr
    PVE_u[[i]] <- Reduce(`+`, ckernel[[i]]$PVE_u)
    if(Reduce(`+`, ckernel[[i]]$PVE_u )<0.01){
      ckernel[[i]] <- NULL
   #   cmat_add[[i]] <- NULL
      cmat_kernel[[i]] <- NULL
      wr[[i]] <- NULL
    }
  }

  lst_u <- sort(unlist(PVE_u), index.return=TRUE, decreasing=TRUE)
  len <- min(length(PVE_u), 50) # set num of 50
  random_index <- lapply(lst_u, `[`, lst_u$x %in% head(lst_u$x,len) ) # select 50 top genes

  if(is.null(cmat_kernel)){
    cmat_kernel_out = cmat_kernel[-which(sapply(cmat_kernel, is.null))]
  }else{
    cmat_kernel_out = cmat_kernel
  }

  cmat_kernel_out = cmat_kernel_out[random_index$ix]

    if(is.null(wr)){
      wr = wr[-which(sapply(wr, is.null))]
  }else{
      wr = wr
  }

  wr = wr[random_index$ix]
  PVE_u = PVE_u[random_index$ix]
#  cmat_kernel_out <- cmat_kernel
  cmat_add_dim = cmat_add_dim[random_index$ix]
 # cmat_add = cmat_add[random_index$ix]
  if(length(cmat_add)>1){
  #  cmat_add_out <- cbind( do.call(cbind, cmat_add), do.call("cbind",cmat_kernel_out) )
    cmat_add_out <- do.call(cbind, cmat_add)
  }else{
  #  cmat_add_out <- cbind(as.matrix(cmat_add), do.call("cbind",cmat_kernel_out) )
    cmat_add_out <- cbind(as.matrix(cmat_add) )
  }

#  cmat_add_out <- cbind(data.sel.geno, data.sel.exp, data.sel.meth)

  return(list(cmat_kernel = cmat_kernel_out, mat_add = cmat_add_out, wr = wr, cmat_add_list = cmat_add_dim,
              PVE_u = PVE_u, gene_sel_index = random_index, cmat_add = cmat_add))
}


