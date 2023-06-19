#################################################
#                                               #
#                Define functions               #
#                                               #
#################################################

####################### Function 0

# ## Function 1
# predictive model

#--------------------------------------------------------------------#

#                            Prediction
#--------------------------------------------------------------------#


vb_predictive <- function(fit, y = NULL,epsilon = 1e-3, B = "No", verbose = FALSE, maf_beta, maf_u){

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

  geno_fix <- fit$genotype0
  # all data - fixed term

  geno_random <- lapply(1:M, function(m) attr(fit,"UD")[[m]] ) # all data - random term

  geno_background <- fit$gamma_lambda # all data - background term

  X <- fit$X  # data for training - fixed term

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

  # calculate and store the PVE_u
  PVE_u <- sapply(1:M, function(m) var_u[[m]]/(post_sigma2e + var_beta + sum(var_u) + post_sigma2b) )

  #PVE_u <- sapply(1:M, function(j) var(TEST$m_u[,j])/var(y) ) # calculate and store the PVE_u

  # sort the PVE_u by the pi_u (decreasing) (add 0 to be first for further calculation)

  sort_PVE_u <- append(0, sapply(1:M, function(m) sum(PVE_u[lst_u$ix][1:m]) ))

  # index of change
  if(length(sort_PVE_u)>1){

    delta_sort_PVE_u <-  sapply(1:M, function(m) delta_sort_PVE_u[m] <- sort_PVE_u[m + 1] - sort_PVE_u[m] )

    index_delta_PVE_u <- which( abs(delta_sort_PVE_u) > epsilon )

    len <- min(length(index_delta_PVE_u), ceiling(0.15*M))

  }else if(length(sort_PVE_u)==1){

    delta_sort_PVE_u <- sort_PVE_u

    index_delta_PVE_u <- 1

    len <- 1
  }

  # finial index for random term # select top regions based on index_delta_PVE_u
  # random_index <- lapply(lst_u, `[`, lst_u$x %in% head(lst_u$x,(index_delta_PVE_u[1] - 1)) )

  random_index <- lapply(lst_u, `[`, lst_u$x %in% head(lst_u$x,len) )

  #------------------------------------------------------------------------------------------------------------
  # (3) calculate the PGE - the proportion of genetic variance explained by the sparse effects terms

  # PGE <- sum(var_beta[fix_index$ix]) / (sum(var_beta[fix_index$ix]) + sum(var_u[random_index$ix]))

  #------------------------------------------------------------------------------------------------------------
  # (4) fixed effects - calculate the predictive dist - fixed term mean

  # data for testing - fixed term
  # X0 <- geno_fix[index_test, index_fixed_marker]
  X0 <- geno_fix[index_test,]

  # Predictive mean
  # all
  # beta_pred0 <- lapply(1:Px, function(j) as.matrix(X0[,j]) %*% w_beta[j] %*% m_beta[j] )

  beta_pred0 <- lapply(1:Px, function(j) as.matrix(X0[,j]) %*% m_beta[j] )

  # extract the selected markers for pred
  beta_pred00 <- beta_pred0[fix_index]

  if (length(fix_index)!=0){
    beta_pred <- Reduce(`+`, beta_pred00) + rep(0, length(index_test))
  }else{
    beta_pred <- as.vector(rep(0, length(index_test)) )
  }

  # (4) random effects - calculate the predictive dist - random term mean

  # data for testing - random term
  # data for testing - random term
  Z0 <- lapply(1:M, function(m) geno_random[[m]][-index_test,] )

  Z1 <- lapply(1:M, function(m) geno_random[[m]][index_test,] )

  # Predictive mean
  # all
  # u_pred0 <- lapply(1:M, function(m) Z0[[m]] %*% (drop(w_u[[m]]) * diag(1, P[[m]])) %*% m_u[[m]] )

  # u0
  u_0 <- lapply(1:M, function(m) Z0[[m]] %*% m_u[[m]] )

  # beta0 for u - prediction step 1

  u0_be <- lapply(1:M, function(m) ginv(crossprod(Z0[[m]])) %*% t(Z0[[m]]) %*% u_0[[m]] )

  # u hat for u - prediction step 2

  u_pred0 <- lapply(1:M, function(m) Z1[[m]] %*% u0_be[[m]] )

  # extract the selected regions for pred
  u_pred00 <- u_pred0[random_index$ix]

  if (length(random_index$ix)!=0){
    u_pred <- Reduce(`+`, u_pred00) + rep(0, length(index_test))
  }else{
    u_pred <- as.vector( rep(0,length(index_test)))
  }

  # (4) background  - calculate the predictive dist - background term mean

  if (B == "Yes"){
    background_pred <- geno_background[index_test,] %*% m_b
  }else{
    background_pred <- as.vector( rep(0,length(index_test)))
  }

  # Predictive variance
  # s_pred <- sqrt(1/model$lambda + diag(X_test %*% model$S %*% t(X_test)))
  #------------------------------------------------------------------------------------------------------------
  # (5) calculate the correlation - accuracy

  pred_value = scale(beta_pred + u_pred + background_pred, T, T)

  m=center
  s=scale #sqrt(scale)  # (scale/sqrt(var(pred_value)))

  pred = ( drop(s) * scale(pred_value) ) + m

#  pred = pred_value

  n_test <- length(index_test)

  y_hat <- mvrnorm(n = 1, mu = as.numeric(beta_pred + u_pred + background_pred) , Sigma = drop(abs(1 - post_sigma2e)/100)*diag(1, n_test) )

  #------------------------------------------------------------------------------------------------------------
  # (6) calculate the correlation - accuracy


  #------------------------------------------------------------------------------------------------------------
  # (6) print the summary info


  if (verbose && is.null(y)==FALSE){

    cat("-------------------------------------------------------------------","\n")
    cat("Variable selection method       :",'Spike and slab',"\n")
    # cat("Sample size for prediction      :",dim(Z0[[1]])[1],"\n")
    cat("No. genetic marker(s)           :",Px,"\n")
    cat("No. genomics region(s)          :",M,"\n")
    cat("No. of top genetic marker(s)    :",length(fix_index$ix),"\n")
    cat("No. of top genomics region(s)   :",length(random_index$ix),"\n")

    cat("\n\n")
    cat("---> Top genetic markers(s):\n")
    for (m in 1:length(fix_index$ix)) {
      cat("Genetic marker(s) with large effects:",fix_index$ix[m],"\n")
    }

    cat("\n\n")
    cat("---> Top region(s):\n")
    for (m in 1:length(random_index$ix)) {
      cat("Region:",random_index$ix[m],"\n")
    }

    #cat("\n\n")
    #cat("---> PGE:",PGE,"\n")

    #cat("\n\n")
    #cat("---> Pearson correlation:",accuracy,"\n")

    #cat("-------------------------------------------------------------------","\n")
  }
  #pred =  beta_pred + u_pred + background_pred
  if(is.null(y)==FALSE){
    return(list(beta_pred = beta_pred, u_pred = u_pred, background_pred = background_pred, pred =  pred,
                PVE_beta = PVE_beta, PVE_u = PVE_u, sort_PVE_u = sort_PVE_u,
                delta_sort_PVE_beta = delta_sort_PVE_beta, delta_sort_PVE_u = delta_sort_PVE_u, y_hat = y_hat,
                fix_index = fix_index, random_index = random_index$ix)) # , s_pred = s_pred
  }else if (is.null(y)==TRUE){
    return(list(beta_pred = beta_pred, u_pred = u_pred, background_pred =background_pred, pred =  pred,
                PVE_beta = PVE_beta, PVE_u = PVE_u, sort_PVE_u = sort_PVE_u,
                delta_sort_PVE_beta = delta_sort_PVE_beta, delta_sort_PVE_u = delta_sort_PVE_u, y_hat = y_hat,
                fix_index = fix_index, random_index = random_index$ix)) # , s_pred = s_pred
  }

}



#add_add_gene <- function(geno, index_test, yc){
#  screening_fix <- bess(geno[-index_test,], yc, family = "gaussian", normalize = TRUE)
#  screening_fix <- bess.one(tot_snp_only0[-index_test,], yc, s = 25, family = "gaussian", normalize = TRUE)
# obtain the top markers
#  subset_fix <- screening_fix$bestmodel$model$xbest
# obtain the index of top markers
#  index_fixed_marker <- attr(subset_fix,"dimnames")[[2]]
#  index_fixed_marker2 <- as.numeric(gsub("[^0-9.]", "",  index_fixed_marker))
#  geno_add <- geno[,index_fixed_marker2]
#  return(geno_add)
#}


## Function 6
# Fit MultiBLUPVariational Bayes model
#Yalu comment out vb_fit_STEP2 <- function(y = y_na, mat_add, kernel_ran,  genoback, maf.filter.snp = 0.01, max_iter = 25000,  weight.type = NULL,
vb_fit_STEP2 <- function(y, mat_add, kernel_ran,  genoback, maf.filter.snp = 0.01, max_iter = 25000,  weight.type = NULL,
                          epsilon_conv = 1e-4, Bsigma2beta = 1, theta_beta = 0.1, theta_u = 0.1, verbose = TRUE){

  #################################################
  #                                               #
  #               0.  Initial Input               #
  #                                               #
  #################################################
  if (verbose ) { cat("[1] Initializing the algorithm:","\n")}
  telliter <- 100

  #--------------------- A. import data --------------------------
  # options(digits=40)

  # Number of observations
  N <- NROW(y)

  # index of the NA value in Y (index of testing data)
  index_test <- as.vector(which(is.na(y)))

  # genotype_sd <-  standardiseGenotypes_sd(genotype0)  # get constant C which is variance of each SNPs

  # 2. matrix normalization

  # genotype <- scale(genotype, center= T, scale=genotype_sd)

  # genotype0 <- scale(genotype00, center= T, scale=T)

  #-------------------------- B. screening for random term and background --------------------------
  # Best subset selection
  #  screening_random <- bess.one(genotype0[-index_test,], y[-index_test], s = nr, family = "gaussian", normalize = TRUE)

  # fix_screening_fix <- bess(genotype_fix, y, family = "gaussian", normalize = TRUE)
  # obtain the top markers
  # subset_random <- screening_random$bestmodel$model$xbest0
  # obtain the index of top markers
  # index_random_marker <- attr(subset_random,"dimnames")[[2]]

  # genotype <- genotype0[,index_random_marker]

  #  genotype <- genotype0

  #genotype <- lapply(1:M, function(m)  genotype[[m]][ , colSums(is.na(genotype[[m]])) == 0] )
  M <- length(kernel_ran)

  Cov_List <- kernel_ran

  # genotype0 <- lapply(1:M0, function(m)  replace(genotype0[[m]],is.na(genotype0[[m]]), 0  ) )
 # genotype0 <- lapply(1:M0, function(m)  genotype0[[m]][ , colSums(is.na(genotype0[[m]])) == 0] )

  # y for training without NA
   yc <- scale(y, center = T, scale = T)[-index_test]

  center <- attr(scale(y, center = T, scale = T),"scaled:center")
  scale <- attr(scale(y, center = T, scale = T),"scaled:scale")

  #--------------------------- select the gene start------------------------------
  # gene_file=genotype0; y=y;

  # ind_gene <- select_gene(gene_file = genotype0, y = y, weight.type)

  # if(is.null(weight.type) == FALSE){
  #   w_g <- lapply(1:M0, function(m) weight_function(genotype0[[m]], N = N, type = weight.type ) )

  # }

  # gene_file = genotype0; y = y; N = N;w = w_g;

  # ind_gene <- select_genesis(gene_file = genotype0, y = y, N = N, w = w_g)

  # var_ind <- sort(ind_gene, index.return=TRUE, decreasing=TRUE)

  # sel_gene_index <- var_ind$ix[1:2] # s 3

  # sel_gene_index <- which(ind_gene > 0) # s 3

 ######

  # genotype <- genotype0[sel_gene_index]

 # genotype <- genotype0

  #--------------------------- select the gene end------------------------------

  #  genotype <- genotype0

#  M <- length(genotype)

#  regions_mat <- lapply(1:M, function(m) as.matrix(scale((genotype[[m]]), center= T, scale= F)) )

  # y for training without NA
  # yc <- scale(y, center = F, scale = F)[-index_test]
  # Z (X)
  #ZList <- lapply(1:M, function(m) regions_mat$mat_int[[m]] %*% UListD[[m]] )
  # combine all snps from all genes
  tot_rare_common_var0 <- do.call("cbind", genoback)

  M0 <- length(genoback)

  # center the big matrix with all snps
  #tot_rare_common_var0 <- as.matrix(scale(tot_rare_common_var0, center= F, scale= F))
  tot_rare_common_var0 <- replace(tot_rare_common_var0,is.na(tot_rare_common_var0), 0)
  tot_rare_common_var1 <- as.matrix(scale(tot_rare_common_var0, center= T, scale= F))

  # DIM <- P  # dim of SNP chunk
  # D <- round(NCOL(genotype)/M)

  # P num of SNPs in each region
  # P <- D

#-------------------------- B. screening for fixed term --------------------------

 #  index_var_snp <- lapply(1:M0, function(m) which(getmafmat(genotype0[[m]]) > maf.filter.snp) )

 # index_var_snp <- lapply(1:M0, function(m) which(getmafmat(genoback[[m]]) > maf.filter.snp) )

 # obtain the common snps by using index

# var_snp <- lapply(1:M0, function(m) genoback[[m]][, index_var_snp[[m]] ] )
# var_snp <- lapply(1:M0, function(m) regions_mat[[m]][, index_var_snp[[m]] ] )

# combine all snps from all genes
#  tot_snp_only0 <- do.call("cbind", var_snp)

  attr(mat_add, "dimnames") <- NULL

 # tot_snp_only0 <- as.matrix(mat_add)

  tot_snp_only0 <- as.matrix(scale(mat_add, center= F, scale= F))

 # tot_snp_only0 <- readRDS(file = "/scale_wlg_persistent/filesets/project/nesi00464/test05062019/snp_condidate.rds")

  #X <- scale(genotype[-index_test,], center= T, scale=standardiseGenotypes_sd(genotype[-index_test,]))

  # screening ---------------------------------------
  #source("G:/My Drive/UoA/PhD_software/VIP/simulation3/gene30/test/CT.R")
  Np_total <- dim(tot_snp_only0)[2]
  Np <- min(N, Np_total, 10)
  # Best subset selection
#  screening_fix <- bess.one(tot_snp_only0[-index_test, keep.var], yc, s = Np, family = "gaussian", normalize = TRUE)
  screening_fix <- bess.one(tot_snp_only0[-index_test, ], yc, s = Np, family = "gaussian", normalize = TRUE)
  # screening_fix <- bess(tot_snp_only0[-index_test,], yc, family = "gaussian", normalize = TRUE)
  # obtain the top markers
  subset_fix <- screening_fix$bestmodel$model$xbest
  # obtain the index of top markers
  index_fixed_marker <- attr(subset_fix,"dimnames")[[2]]
  index_fixed_marker2 <- as.numeric(gsub("[^0-9.]", "",  index_fixed_marker))
  # get the index of common snps maf > 0. 01
  # index_var_snp <- lapply(1:M, function(m) which(maf(t(genotype[[m]])) > maf.filter.snp) )

  # var_maf <- abs(maf(t(Z)))

  # tot_snp_only <- readRDS(file = "/scale_wlg_persistent/filesets/project/nesi00464/test05062019/snp_condidate.rds")

  # X <- subset_fix
  tot_snp_only <- tot_snp_only0[,index_fixed_marker2]
 # tot_snp_only <- tot_snp_only0[,keep.var]
 # tot_rare_common_noC <- tot_rare_common_var0[, index_fixed_marker2]
 # tot_rare_common_C <- tot_rare_common_var1[, index_fixed_marker2]

# tot_snp_only <- mat_add

#  tot_snp_only <- tot_rare_common_C

  X <- tot_snp_only[-index_test,]

  #  X <- as.matrix(scale(X, center= T, scale= F))

  Px <- round(NCOL(X))  # P num of SNPs in each region

  ##-------------------------- add gene_add to genotype ----------------------------------

  # geno=tot_snp_only0; index_test=index_test; yc=yc;

  # gene_add <- add_add_gene(geno=tot_snp_only0, index_test=index_test, yc=yc)

  #  genotype[[length(genotype)+1]] <- gene_add

  #  M <- length(genotype)

  #  regions_mat <- lapply(1:M, function(m) as.matrix(scale(genotype[[m]], center= T, scale= F)) )

  ##------------------------------------------------
  ## get the weight for each genetic region

 # if(is.null(weight.type) == FALSE){

#    w0 <- lapply(1:M, function(m) weight_function(genotype[[m]], N = N, type = weight.type ) )

#    if(weight.type=="uw"){
#      w  <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))) )
#      pve <- 0.8
#    }else if(weight.type=="beta"){
#      w  <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))) )
#      pve <- 0.8
#    }else if(weight.type=="wss"){
#      w  <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))) )
#      pve <- 0.8
#    }

    # data=genotype[[2]]; N = N; type = weight.type;

#  }

#  gc(verbose = FALSE)

  #-------------------------- C. reconstruct the covarience matrix for random term--------------------------

  # 1.  Z original SNP matrix
  # ZList0 <- lapply(1:M, function(m)  dropNA(regions_mat[[m]])  )

  #  ZList0 <- lapply(1:M, function(m)  regions_mat[[m]] )

#  NP <- lapply(1:M, function(m)  dim(regions_mat[[m]])[2] )
#  ZList0 <- lapply(1:M, function(m)  scale(regions_mat[[m]], F, F) )

  # A_ZList_sd <- lapply(1:M, function(m) standardiseGenotypes_sd(ZList00[[m]]))  # get constant C which is variance of each SNPs

  # 2. matrix normalization

  # ZList0 <- lapply(1:M, function(m) scale(ZList00[[m]], center= T, scale=A_ZList_sd[[m]]) )

  #ZList0 <- ZList00

  # ZList0 <- ZList00 # for simulation data

  # 3. select the kernel function - Lin, RBF, IBS, POLY2
  ## weight function - uw, wss, beta
  ## calculate kernel matrix

  #if(kernel=="Lin" && is.null(weight.type) == TRUE){
  #  Cov_List <- lapply(1:M, function(m) tcrossprod(ZList0[[m]]) )  # Linear effects for SNPs
  #}else if(kernel=="RBF" && is.null(weight.type) == TRUE){
  #  Cov_List  <- lapply(1:M, function(m) (kernelMatrix(rbfdot(sigma = 0.3), ZList0[[m]])@.Data) ) # RBF effects for SNPs
  #}else if(kernel=="POLY2" && is.null(weight.type) == TRUE){
  #  Cov_List <- lapply(1:M, function(m) (tcrossprod(ZList0[[m]]) + 1)**2 ) # Poly 2 effects for SNPs
  #}else if(kernel=="Lin" && is.null(weight.type) == FALSE ){
  #  Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]]))/NP[[m]] )
  #  Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]])) )
  #}else if(kernel=="quadratic" && is.null(weight.type) == FALSE ){
  #  Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2)/NP[[m]] )
  #  Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2) )
  #}else{
  #  stop ("currently only support Lin, RBF, POLY2 and IBS kernel and uw, beta, wss weight.type")
  #}

  #------------------------------- pc selection -----------------------------

  # reconstruct the covariance matrix using eigenvalues and eigen vectors - random effects

  # Cov_List2 <- lapply(1:M, function(m) cov2cor(Cov_List[[m]]) )

  # EigenList <- lapply(1:M, function(m) eigen(Cov_List[[m]]) ) # get eigen values and vectors of the covariance

  #  EigenvalueList0 <- lapply(1:M, function(m) ifelse(EigenList[[m]]$values < 0, 0, EigenList[[m]]$values) ) # set negative eigen values to zero

  #DIM_N <- lapply(1:M, function(m) first(which(EigenList[[m]]$values < 0 )) ) # set negative eigen values to zero

  # vectorlist0 <- lapply(1:M, function(m) EigenList[[m]]$vectors[,1:DIM_N[[m]]] ) # extract the eigen vectors with top eigen values
  #  vectorlist0 <- lapply(1:M, function(m) EigenList[[m]]$vectors ) # extract the eigen vectors with top eigen values

  # ZList00 <- lapply(1:M, function(m) vectorlist0[[m]] %*% diag(EigenvalueList0[[m]][1:DIM_N[[m]]]) %*% t(vectorlist0[[m]]) ) # reconstruct X - (beta = U * D^(1/2))
  #  ZList00 <- lapply(1:M, function(m) vectorlist0[[m]] %*% diag(EigenvalueList0[[m]]) %*% t(vectorlist0[[m]]) ) # reconstruct X - (beta = U * D^(1/2))

  for(i in 1:M){
    attr(Cov_List[[i]], "dimnames") <- NULL
  } # remove attr

  # a. Proportion of Variance
  #  temp <- lapply(1:M, function(m) princomp(ZList00[[m]],cor = TRUE, scores=TRUE) )
  temp <- lapply(1:M, function(m) prcomp(Cov_List[[m]]) )
  zv <- lapply(1:M, function(m) temp[[m]]$sdev )
  pv <- lapply(1:M, function(m) zv[[m]]^2 )
  cs <- lapply(1:M, function(m) pv[[m]]/sum(pv[[m]]) )
  # b. Cumulative Proportion
  cm <- lapply(1:M, function(m) cumsum(cs[[m]]) )
  # c. select number of pc explaining 80% of variance
#  DIM <- lapply(1:M, function(m) first(which(cm[[m]]>pve)) )

  # DIM <- lapply(1:M, function(m) min(first(which(cm[[m]]>pve)), 15) )

  # t.wss3 <- Cov_List <- lapply(1:M, function(m) tr(Cov_List[[m]]) )

  #  unlist(t.uw)
  #  unlist(t.beta)
  #  unlist(t.wss)
  #  unlist(t.wss2)
  #  unlist(t.wss3)

  # 4. reconstruct the covariance matrix using eigenvalues and eigen vectors - random effects

  # Cov_List2 <- lapply(1:M, function(m) cov2cor(Cov_List[[m]]) )

  Cov_List <- lapply(1:M, function(m) (Cov_List[[m]])/(tr(Cov_List[[m]])*1) )

  EigenList <- lapply(1:M, function(m) eigen(Cov_List[[m]]) ) # get eigen values and vectors of the covariance

  EigenvalueList0 <- lapply(1:M, function(m) ifelse(EigenList[[m]]$values < 0, 0, EigenList[[m]]$values) ) # set negative eigen values to zero

  #  EigenvalueList0 <- lapply(1:M, function(m) scale(EigenvalueList0[[m]], F, scale = T) ) # set negative eigen values to zero

  DIM <- lapply(1:M, function(m) max(min(length(which(EigenvalueList0[[m]] >1e-30)), 10), 1) )

  indexlist <- lapply(1:M, function(m) sort(EigenvalueList0[[m]], index.return=TRUE, decreasing=TRUE) ) # get the eigen values sorted index

  vectorlist <- lapply(1:M, function(m) EigenList[[m]]$vectors[,indexlist[[m]]$ix[1:DIM[[m]]]] ) # extract the eigen vectors with top eigen values

  EigenvalueList <- lapply(1:M, function(m) indexlist[[m]]$x[1:DIM[[m]]] ) # get the top eigen values index (Top P)

  # U * (D)^1/2

  # value.uw <- lapply(1:M, function(m) sum(EigenvalueList[[m]]))

  # unlist(value.beta)

  UListD <- lapply(1:M, function(m) (vectorlist[[m]] %*% diag(sqrt(EigenvalueList[[m]][ 1:DIM[[m]] ]), DIM[[m]]))/10 ) # reconstruct X - (beta = U * D^(1/2))

  # reconstructed Z matrix

  Z <- lapply(1:M, function(m) UListD[[m]][-index_test,] )

  tUDUList <- lapply(1:M, function(m) crossprod(UListD[[m]][-index_test,]) ) # t(U)*D*U

  # reconstructed covariance

  tZZList <- tUDUList

  #-----------------------------------------------------------------------------------------------------------------
  # 5. reconstruct the covariance matrix using eigenvalues and eigen vectors - background

  # Cov_b <- tcrossprod(tot_snp_only) # only snps

  # tot_var <- dropNA(tot_rare_common_var) # dropNA all variants

  NP0 <- dim(tot_rare_common_var0)[2]

  # w_g0 <- weight_function(tot_rare_common_var0, N = N, type = weight.type )

  #  w_g  <- diag(w_g0)

  tot_var <- tot_rare_common_var0 # dropNA all variants

  #  Cov_b <- (tot_var %*% w_g %*% t(tot_var)) # all snps

  if(is.null(weight.type) == FALSE){

    weight.type="uw"

    w_g0 <- weight_function(tot_rare_common_var0, N = N, type = weight.type )

    if(weight.type=="uw"){
      w_g  <-  diag(as.numeric(scale(w_g0, F, scale = F)))
    }else{
      w_g  <-  diag(as.numeric(scale(w_g0, F, scale = T)))
    }

    # data=genotype[[2]]; N = N; type = weight.type;

  }

  Cov_b <- (tot_var %*% w_g %*% t(tot_var))/NP0 # all snps

  #  Cov_b <- (tot_var %*% t(tot_var))/NP0 # all snps

  Eigen_b <- eigen(Cov_b) # get eigen values and vectors of the covariance

  Eigenvalue_b0 <- ifelse(Eigen_b$values < 0, 0, Eigen_b$values) # set negative eigen values to zero

  index_b <-  sort(Eigenvalue_b0, index.return=TRUE, decreasing=TRUE)  # get the eigen values sorted index

  # DIM0 <- length(which(Eigenvalue_b0!=0)) # remove eigenvalue 0

  # a. Proportion of Variance for g0 ----------------------------
  #  temp <- lapply(1:M, function(m) princomp(ZList00[[m]],cor = TRUE, scores=TRUE) )
  temp0 <- prcomp(Cov_b)
  zv0 <-  temp0$sdev
  pv0 <-  zv0^2
  cs0 <-  pv0/sum(pv0)
  # b. Cumulative Proportion
  cm0 <- cumsum(cs0)
  # c. select number of pc explaining 80% of variance
  DIM0 <- max(first(which(cm0 >0.9)),2)

  # DIM0 <- max(min(length(which(Eigenvalue_b0 >1e-30)), 15), 1)

  vector_b <- Eigen_b$vectors[,index_b$ix[1:DIM0]]  # extract the eigen vectors with top eigen values

  Eigenvalue_b <- index_b$x[1:DIM0] # get the top eigen values index (Top P)

  # Gamma * (Lambda)^1/2

  gamma_lambda <- vector_b %*% diag(sqrt(Eigenvalue_b[1:DIM0])) # reconstruct X - (beta = U * D^(1/2))

  # reconstructed g matrix

  gamma_lambda_b <- gamma_lambda[-index_test,]  # gamma*(lambda)^(1/2)

  tglg <- crossprod(gamma_lambda_b)  # t(gamma)*lambda*gamma

  # reconstructed covariance - background

  # tZZList <- tUDUList

  gc(verbose = FALSE)

  #--------------------------- D. initial the parameters for iterations ----------------------------------

  #s_beta

  #E_theta <- lapply(1:M, function(m) A_theta_u0[[m]]/(A_theta_u0[[m]] + B_theta_u0[[m]]) )

  # fixed theta and sigma of beta for fixed term
  #### fix ------------------------
  A_sigma2beta0 <- 0.1
  B_sigma2beta0 <- 0.1

  A_sigma2beta <- 0.1
  B_sigma2beta <- 0.1

  A_theta_beta0 <- 0.1
  B_theta_beta0 <- 0.1

  E_theta_beta <- rep(0.1, Px)

  B_sigma2beta  <- 0.1 # fix

  #m_beta

  m_beta <- rep(1/(Px*N), Px)

  E_w_beta <- rep(0.1, Px)

  XGammaBeta <- X %*% diag(E_w_beta) %*% m_beta

  logit_gamma_beta <- rep((0.05/0.95), Px)

  E_w_beta <- sapply(1:Px, function(j) inv.logit(logit_gamma_beta[j]) )

  omega <- tcrossprod(E_w_beta) + hadamard.prod(diag(E_w_beta), (diag(1,Px)-diag(E_w_beta)))

  #B_sigma2beta <- rep(0.1, Px)

  #### random effects ------------------------

  A_theta_u0 <- 0.1
  B_theta_u0 <- 0.1

  A_sigma2u0 <- rep(0.1, M)
  B_sigma2u0 <- rep(0.1, M)

  A_sigma2u <- rep(A_sigma2u0 + N/2, M)
  B_sigma2u <- rep(0.1, M)

  m_u <- lapply(1:M, function(m) rep(0.1 , DIM[[m]]) )

  logit_gamma_u <- rep(0.05/0.95, M)

  E_w_u <- lapply(1:M, function(m) inv.logit(logit_gamma_u[[m]]) )

  E_theta_u <- rep(0.1, M)

  ZGammaU0 <- lapply(1:M, function(m) Z[[m]] %*% (drop(E_w_u[[m]]) * diag(1, DIM[[m]])) %*% m_u[[m]]
  )

  ZGammaU <- do.call("cbind", ZGammaU0)

  # E_theta_u <- rep(E_theta_u00,M) ## random

  ####  background ------------------------
  m_b <-  rep(0.1, DIM0)

  A_sigma2b <- 0.1
  B_sigma2b <- 0.1

  A_sigma2b0 <- 0.1
  B_sigma2b0 <- 0.1

  g <- gamma_lambda_b %*% m_b

  ####  error ------------------------

  A_sigma2e0 <- 0.1
  B_sigma2e0 <- 0.1

  A_sigma2e <- A_sigma2e0 + N/2
  B_sigma2e <- 0.1

  #----------------------------------------------------------------------
  #----------------------------------------------------------------------
  #----------------------------------------------------------------------

  # Store the lower bounds L
  L <- rep(-Inf, max_iter)
  if (verbose ) { cat("[1] Done!","\n")}
  if (verbose ) { cat("[2] Starting variational inference:","\n")}

  # Iterate to find optimal parameters
  for (i in 2:max_iter) {

    #################################################
    #                                               #
    #      1.  Initial Paramters and Iteration      #
    #                                               #
    #################################################

    #----------------------# random effects #----------------------#

    B   <- sapply(1:M, function(m) yc - XGammaBeta - rowSums(ZGammaU[,-m, drop=F]) - g
    )

    # Um

    s_u <- lapply(1:M, function(m) solve( drop(((A_sigma2e0 + N/2)/B_sigma2e)) * hadamard.prod(crossprod(Z[[m]]), (drop(E_w_u[[m]]) * diag(1, DIM[[m]]))) +
                                            drop(((A_sigma2u0[[m]] + DIM[[m]]/2))/B_sigma2u[[m]]) * diag(1, DIM[[m]])
    )
    )

    m_u <- lapply(1:M, function(m) (drop((A_sigma2e0 + N/2)/B_sigma2e) * s_u[[m]]) %*% ((drop(E_w_u[[m]])) * diag(1, DIM[[m]])) %*% t(Z[[m]]) %*% B[,m, drop=F]
    )

    # Gamma w_u

    logit_gamma_u <- lapply(1:M, function(m) logit(E_theta_u[m]) + drop((A_sigma2e0 + N/2)/B_sigma2e) * (
      t(B[,m, drop=F]) %*% Z[[m]] %*% m_u[[m]] - (1/2) * (
        t(m_u[[m]]) %*% t(Z[[m]]) %*% Z[[m]] %*% m_u[[m]] +
          tr(crossprod(Z[[m]]) %*% s_u[[m]])
      )
    )
    )

    E_w_u <- lapply(1:M, function(m) inv.logit(logit_gamma_u[[m]])
    )

    ZGammaU <- sapply(1:M, function(m) Z[[m]] %*% (drop(E_w_u[[m]]) * diag(1, DIM[[m]])) %*% m_u[[m]]
    )

    # sigma B of U_m

    B_sigma2u <- lapply(1:M, function(m) B_sigma2u0[[m]] + (1/2) * ( crossprod(m_u[[m]]) + tr(s_u[[m]])
    )
    )

    #  B_sigma2u <- lapply(1:M, function(m) 0.1 )

    # theta for random effects

    A_theta_u  <- sapply(1:M, function(m) E_w_u[[m]] + A_theta_u0
    )

    B_theta_u  <- sapply(1:M, function(m) 1 - E_w_u[[m]] + B_theta_u0
    )

    E_theta_u  <- sapply(1:M, function(m) A_theta_u[m]/(A_theta_u[m] + B_theta_u[m])
    )
    #  E_theta_u  <- rep(theta_u,M)

    #----------------------# fixed effects #----------------------#

    A      <- yc - rowSums(ZGammaU[,, drop=F]) - g

    # Beta

    s_beta <- solve( drop(((A_sigma2e0 + N/2)/B_sigma2e)) * hadamard.prod(crossprod(X), omega) +
                       drop((A_sigma2beta0 + Px/2)/B_sigma2beta) * diag(1, Px)
    )

    m_beta <-  drop(((A_sigma2e0 + N/2)/B_sigma2e)) * s_beta %*% diag(E_w_beta) %*% t(X) %*% A


    # sigma2 for fixed effects

    A_sigma2beta <- (1/2) * N + A_sigma2beta0

    mtm_beta <-  crossprod(m_beta)

    tr_s_beta <- tr(s_beta)

    B_sigma2beta <- (1/2) * (mtm_beta + tr_s_beta) + B_sigma2beta0

    # B_sigma2beta <- B_sigma2beta0

    # Gamma w_beta

    # use X_E_w_beta[,-j] = X[,-j] %*% diag(E_w_beta[-j]) to replace (for computation purpose)
    X_E_w_beta <- X %*% diag(E_w_beta)

    logit_gamma_beta <- sapply(1:Px, function(j) logit(E_theta_beta[j]) -
                                 (1/2) * drop((A_sigma2e0 + N/2)/B_sigma2e) * crossprod(X[,j]) %*% (crossprod(m_beta[j]) + s_beta[j,j]) +
                                 drop(((A_sigma2e0 + N/2)/B_sigma2e)) * t(X[,j]) %*% ( A * m_beta[j] -
                                                                                         X_E_w_beta[,-j] %*% (m_beta[-j] * m_beta[j] + s_beta[-j,j])
                                 )
    )

    E_w_beta <- sapply(1:Px, function(j) inv.logit(logit_gamma_beta[j])
    )

    # theta for fixed effects

    A_theta_beta <- sapply(1:Px, function(j) E_w_beta[j] + A_theta_beta0
    )

    B_theta_beta  <- sapply(1:Px, function(j) 1 - E_w_beta[j] + B_theta_beta0
    )

    E_theta_beta  <- sapply(1:Px, function(j) A_theta_beta[j]/(A_theta_beta[j] + B_theta_beta[j])
    )

    #  E_theta_beta <- rep(theta_beta, Px)

    # Omega
    omega  <- tcrossprod(E_w_beta) + hadamard.prod(diag(E_w_beta), (diag(1,Px)-diag(E_w_beta)))

    XGammaBeta <- X %*% diag(E_w_beta) %*% m_beta

    #----------------------# background #----------------------#

    B0      <- yc - XGammaBeta - rowSums(ZGammaU[,, drop=F])

    # Beta

    s_b <- solve( drop(((A_sigma2e0 + N/2)/B_sigma2e)) * tglg +
                    drop((A_sigma2b0 + DIM0/2)/B_sigma2b) * diag(1, DIM0)
    )

    m_b <-  drop(((A_sigma2e0 + N/2)/B_sigma2e)) * s_b %*% t(gamma_lambda_b) %*% B0

    A_sigma2b <- (1/2) * N + A_sigma2b0

    mtm_b <-  crossprod(m_b)

    tr_s_b <- tr(s_b)

    B_sigma2b <- (1/2) * (mtm_b + tr_s_b) + B_sigma2b0

    # background g
    g <- gamma_lambda_b %*% m_b


    #----------------------# error #----------------------#

    # sigma B of error

    #  CtC1  <- 2 * t(yc) %*% ( X %*% diag(E_w_beta) %*% m_beta
    #  )

    CtC1  <- X %*% diag(E_w_beta) %*%  m_beta

    CtC2  <- sapply(1:M, function(m) drop(E_w_u[[m]]) * Z[[m]] %*% m_u[[m]]
    )

    CtC3 <- g
    ####
    CtC4  <- tr( hadamard.prod(crossprod(X), omega) %*% s_beta
    )

    CtC5  <- sapply(1:M, function(m) tr( drop(E_w_u[[m]]) * crossprod(Z[[m]]) %*% s_u[[m]]
    )
    )

    CtC6  <- tr(tglg %*% s_b
    )

    B_sigma2e <- B_sigma2e0 + (1/2) * ( crossprod(yc - CtC1 - rowSums(CtC2) - CtC3 ) +
                                          sum(CtC4) +
                                          sum(CtC5) +
                                          sum(CtC6)
    )


    #---------------- calculating the sigma --------------------------

    ##########################################

    post_sigma2beta <- B_sigma2beta/(A_sigma2beta - 1)

    post_sigma2u <- sapply(1:M, function(m) B_sigma2u[[m]] / (A_sigma2u[m] - 1) )

    post_sigma2b <- B_sigma2b / (A_sigma2b - 1)

    post_sigma2e <- B_sigma2e / (A_sigma2e - 1)

    ##########################################


    #################################################
    #                                               #
    #          2.  Iteration (lower Bound)          #
    #                                               #
    #################################################

    # Compute lower bound

    #--------------------p--------------------

    # lb_py <- -N/2 * log(2 * pi) - 1/2 * (A_sigma2e/B_sigma2e) * ( crossprod(yc - CtC1 - sum(CtC2) ) +
    #                                                                 sum(CtC3) + sum(CtC4))

    lb_py <- -N/2 * log(2 * pi) - 1/2 * (A_sigma2e/B_sigma2e) * ( crossprod(yc  - rowSums(CtC2) ) +
                                                                    sum(CtC3) )

    # fix

    lb_pbeta <- N/2 * log(2 * pi) - 1/2 * (A_sigma2beta/B_sigma2beta) * (crossprod(m_beta) + tr(s_beta))

    lb_psigma2_beta <- A_sigma2beta0 * log(B_sigma2beta0) - lgamma(A_sigma2beta0) -
      (A_sigma2beta0 + 1) * (B_sigma2beta/(A_sigma2beta - 1)) + B_sigma2beta0 * (A_sigma2beta/B_sigma2beta)

    lb_pgamma_beta <- sapply(1:Px, function(j) E_w_beta[j] * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
                               (1 - E_w_beta[j]) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
    )

    lb_ptheta_beta <- sapply(1:Px, function(j) (A_theta_beta0 - 1) * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
                               (B_theta_beta0 - 1) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
    )

    # random

    lb_pu <- sapply(1:M, function(m) - N/2 * log(2 * pi) - N/2 * (A_sigma2u[m]/B_sigma2u[[m]]) * (crossprod(m_u[[m]]) + tr(s_u[[m]]))
    )

    lb_psigma2_u <- sapply(1:M, function(m) A_sigma2u0[m] * log(B_sigma2u0[m]) - lgamma(A_sigma2u0[m]) -
                             (A_sigma2u0[m] + 1) * (log(B_sigma2u[[m]]) - digamma(A_sigma2u[m])) +
                             B_sigma2u0[m] * (A_sigma2u[m]/B_sigma2u[[m]])
    )

    lb_pgamma_u <- sapply(1:M, function(m) E_w_u[[m]] * (digamma(A_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m])) +
                            (1 - E_w_u[[m]]) * (digamma(B_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m]))
    )

    lb_ptheta_u <- sapply(1:M, function(m) (A_theta_u0 - 1) * (digamma(A_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m])) +
                            (B_theta_u0 - 1) * (digamma(B_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m]))
    )

    # background

    lb_pb <- N/2 * log(2 * pi) - 1/2 * (A_sigma2b/B_sigma2b) * (crossprod(m_b) + tr(s_b))

    lb_psigma2_b <- A_sigma2b0 * log(B_sigma2b0) - lgamma(A_sigma2b0) -
      (A_sigma2b0 + 1) * (B_sigma2b/(A_sigma2b - 1)) + B_sigma2b0 * (A_sigma2b/B_sigma2b)

    # error

    lb_psigma2_e <- A_sigma2e0 * log(B_sigma2e0) - lgamma(A_sigma2e0) -
      (A_sigma2e0 + 1) * (log(abs(B_sigma2e)) - digamma(A_sigma2e)) -
      B_sigma2e0 * (A_sigma2e/B_sigma2e)

    #--------------------q--------------------

    # fix

    lb_qsigma2_beta <- - A_sigma2beta - log(B_sigma2beta) * lgamma(A_sigma2beta) +
      (1 + A_sigma2beta) * digamma(A_sigma2beta)

    #  lb_qbeta <- - 1/2 * log(exp(determinant(s_beta)$modulus[1]) + 1e-10) - 1/2 * log(2 * pi)

    lb_qbeta <- - 1/2 * determinant(s_beta)$modulus[1] - 1/2 * log(2 * pi)

    lb_qgamma_beta <- sapply(1:Px, function(j) E_w_beta[j] * log(E_w_beta[j] + 1e-10) +
                               (1 - E_w_beta[j]) * log(1 - E_w_beta[j] + 1e-10)
    )

    lb_qtheta_beta <- sapply(1:Px, function(j) lgamma(A_theta_beta[j] + B_theta_beta[j]) - lgamma(A_theta_beta[j]) - lgamma(B_theta_beta[j]) +
                               (A_theta_beta[j] - 1) * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
                               (B_theta_beta[j] - 1) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
    )

    # random

    lb_qsigma2_u <- sapply(1:M, function(m) - A_sigma2u[m] - (log(B_sigma2u[[m]]) + log(lgamma(A_sigma2u[m]))) + (1 + A_sigma2u[m]) *
                             digamma(A_sigma2u[m]) )

    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(exp(determinant(s_u[[m]])$modulus[1]) + 1e-10) -1/2 * log(2 * pi) )

    lb_qu <- sapply(1:M, function(m) - 1/2 * determinant(s_u[[m]])$modulus[1] -1/2 * log(2 * pi) )
    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(abs(determinant(s_u[[m]])$modulus[1]) + 1e-10) -1/2 * log(2 * pi) )

    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(det(s_u[[m]]) + 1e-10) -1/2 * log(2 * pi) )
    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(ifelse(is.infinite(det(s_u[[m]])), 1e20, det(s_u[[m]])) + 1e-10) -1/2 * log(2 * pi) )

    lb_qgamma_u <- sapply(1:M, function(m) E_w_u[[m]] * log(E_w_u[[m]] + 1e-10) +
                            (1 - E_w_u[[m]]) * log(1 - E_w_u[[m]] + 1e-10)
    )

    lb_qtheta_u <- sapply(1:M, function(m) lgamma(A_theta_u[m] + B_theta_u[m]) - lgamma(A_theta_u[m]) - lgamma(B_theta_u[m]) +
                            (A_theta_u[m] - 1) * (digamma(A_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m])) +
                            (B_theta_u[m] - 1) * (digamma(B_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m]))
    )

    # background

    lb_qsigma2_b <- - A_sigma2b - log(B_sigma2b) * lgamma(A_sigma2b) +
      (1 + A_sigma2b) * digamma(A_sigma2b)

    #  lb_qb <- - 1/2 * log(exp(determinant(s_b)$modulus[1]) + 1e-10) - 1/2 * log(2 * pi)

    lb_qb <- - 1/2 * determinant(s_b)$modulus[1] - 1/2 * log(2 * pi)

    # error
    lb_qsigma2_e <- -A_sigma2e - (log(abs(B_sigma2e)) + log(lgamma(A_sigma2e))) + (1 + A_sigma2e) * digamma(A_sigma2e)

    #---------------------L--------------------

    L[i] <- lb_py +
      lb_pbeta + lb_psigma2_beta + sum(lb_pgamma_beta) + sum(lb_ptheta_beta) +
      sum(lb_pu) + sum(lb_psigma2_u) + sum(lb_pgamma_u) + sum(lb_ptheta_u) +
      lb_pb + lb_psigma2_b +
      lb_psigma2_e  -
      lb_qsigma2_beta - lb_qbeta - sum(lb_qgamma_beta) - sum(lb_qtheta_beta) -
      sum(lb_qu) - sum(lb_qsigma2_u) - sum(lb_qgamma_u) - sum(lb_qtheta_u) -
      lb_qsigma2_b - lb_qb -
      lb_qsigma2_e

    log_delta <- abs(L[i] - L[i - 1])
    # Brute force
    #  L0 <- L
    #  L0 <- L0[is.finite(L0)]

    #  if(length(unique(L0)) == length(L0) && i != max_iter){
    #    log_delta <- abs(L[i] - L[i - 1])
    #  }else if(length(unique(L0)) != length(L0) || i == max_iter){
    #    log_delta <- 1e-10
    #  }

    #################################################
    #                                               #
    #           3. Output Iteration Results         #
    #                                               #
    #################################################

    # Show VB difference
    if (verbose && i %% telliter == 0) { cat("Iteration:\t",i,"\tLB:\t", L[i],"\tDelta:\t",
                                             log_delta,"\n")}

    # Check for convergence
    if (log_delta < epsilon_conv) {

      if (verbose && i %% telliter != 0) { cat("Iteration:\t",i,"\tLB:\t", L[i],"\tDelta:\t",
                                               log_delta,"\n")}

      cat("[2] VB coverged!\n")
      break }
    # Check if VB converged in the given maximum iterations
    if (i == max_iter ) {message("[2] VB did not converge!\n")}

  }

  if (log_delta < epsilon_conv) {
    obj <- structure(list(X = X, Z = Z, genotype0 = tot_snp_only, gamma_lambda = gamma_lambda, B_sigma2beta= B_sigma2beta,
                          index_test = index_test, B_sigma2e = B_sigma2e, B_sigma2u = B_sigma2u, post_sigma2b =post_sigma2b,
                          post_sigma2u = post_sigma2u, post_sigma2beta = post_sigma2beta, post_sigma2e = post_sigma2e,
                          N = N, M = M, Px = Px, DIM = DIM, DIM0 = DIM0, Delta = log_delta, E_theta_u = E_theta_u,
                          E_theta_beta = E_theta_beta, s_beta = s_beta, m_b = m_b, s_b = s_b,
                          m_beta = m_beta,r = M,m_u = m_u,s_u = s_u, y=yc, scale = scale, center = center,
                          E_w_beta = E_w_beta, E_w_u = E_w_u, Px =Px,index_fixed_marker = index_fixed_marker2,
                          L = L[2:i]), UD = UListD, class = "vb")
    return(obj)

  }

}

