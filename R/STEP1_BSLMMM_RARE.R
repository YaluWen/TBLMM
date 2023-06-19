

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

#---------------------------------------------------------------------------------------
Lin0<-function(x) tcrossprod(x)/max(1,ncol(x))

Quad1<-function(x) (tcrossprod(x)+1)^2/max(1,ncol(x))

IBS<-function(x)  1 - as.matrix(dist(x, method='manhattan') * .5 /max(1, ncol(x)) )  ## dist does scaLing in the presence of missing values

Gau <-function(x){p=ncol(x);return(exp(-1*as.matrix(dist(x)^2)/p))}

NNnetwork<-function(x)

{

  L1=Lin0(x)

  L2=diag((tcrossprod(x)/max(1,ncol(x))+1));

  L3=matrix(L2,ncol=1) %*% matrix(L2,nrow=1)

  asin(L1/sqrt(L3))

}

#------------------------------------------------------------#

## Function 2
# get Allele Frequencies
getAlleleFrequencies <- function(snp) {
  if (dim(data.frame(snp))[2] != 1) {
    stop ("getAllelelFrequencies takes a [N x 1] vector of genotypes, but ",
          "provided genotype dimensions are: [", dim(snp)[1], " x ",
          dim(snp)[2], "]")
  }

  pp <- length(which(snp == 0))
  pq2 <- length(which(snp == 1))
  qq <- length(which(snp == 2))

  if (pp + pq2 + qq != length(snp[!is.na(snp)])) {
    # stop ("SNP vector contains alleles not encoded as 0, 1 or 2")
    p <- 0.5
  }else{
    p <- (2*pp +  pq2)/(2*length(snp[!is.na(snp)]))
  }
  return(c(p, 1-p))
}

getmafmat <- function(mat) {

  nsnp = dim(mat)[2]
  var_maf <- sapply(1:nsnp, function(i) min( getmaf(mat[,i]) ) )
  return(var_maf)
}

getmaf <- function(snp) {
  if (dim(data.frame(snp))[2] != 1) {
    stop ("getAllelelFrequencies takes a [N x 1] vector of genotypes, but ",
          "provided genotype dimensions are: [", dim(snp)[1], " x ",
          dim(snp)[2], "]")
  }

  pp <- length(which(snp == 0))
  pq2 <- length(which(snp == 1))
  qq <- length(which(snp == 2))

  if (pp + pq2 + qq != length(snp[!is.na(snp)])) {
   # stop ("SNP vector contains alleles not encoded as 0, 1 or 2")
    p <- 0.5
  }else{
    p <- (2*pp +  pq2)/(2*length(snp[!is.na(snp)]))
  }
  return(min(p, 1-p))
}

## Function 3
# get Genotypes sd
standardiseGenotypes_sd <- function(geno) {
  allele_freq <-  sapply(data.frame(geno),  getAlleleFrequencies)
  var_geno <- sqrt(2*allele_freq[1,]*allele_freq[2,])
  var_geno[var_geno == 0] <- 1
  return (var_geno)
}

## Function 4
# get Genotypes mean
standardiseGenotypes_mean <- function(geno) {
  allele_freq <-  sapply(data.frame(geno),  getAlleleFrequencies)
  var_mean <- 2*allele_freq[1,]
  return (var_mean)
}


## Function 5
# cut into regions

## function 6 weight function
weight_function <- function(data, N, type = "uw"){

  Z = matrix(unlist(data), nrow = N, byrow = F)

  nsnp = dim(Z)[2]

  # var_maf <- abs(maf(t(Z)))

  var_maf <- sapply(1:nsnp, function(i) min( getmaf(Z[,i]) ) )

  #  var_maf <-colMeans(Z, na.rm = TRUE)/2

  if (type == "uw") {

    w = rep(1, length(var_maf))

  }else if (type == "beta") {

    w = dbeta(var_maf, 1, 25)^2

  }else if (type == "wss") {

    w = 1/(var_maf * (1-var_maf))

  }

  w[which(!is.finite(w))] <- 0

  return(w)
}

crossprod.replacena <- function(x, val=0) {
  crossprod(replace(x, is.na(x), val)
  )
}

crossprod.replacena2 <- function(x, y, val=0) {
  crossprod(replace(x, is.na(x), val),
            replace(y, is.na(y), val)
  )
}

tcrossprod.replacena <- function(x, val=0) {
  tcrossprod(replace(x, is.na(x), val)
  )
}

tcrossprod.replacena2 <- function(x, y, val=0) {
  tcrossprod(replace(x, is.na(x), val),
             replace(y, is.na(y), val)
  )
}

## Function 6
# Fit MultiBLUPVariational Bayes model
#' @import gtools
#' @importFrom matrixcalc hadamard.prod
#' @importFrom psych tr
#' @importFrom data.table first
#' @importFrom BeSS bess.one
#' @importFrom stats dist
#' @importFrom stats dbeta
#' @importFrom stats prcomp
#' @importFrom stats var
vb_fit_STEP1 <- function(y, genotype0, max_iter = 1000, weight.type = NULL, maf.filter.snp = 0.01, epsilon_conv = 1e-4,
                              Bsigma2beta = 1, theta_beta = 0.1, theta_u = 0.1, verbose = TRUE, kernel = "Lin"){

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

  #-------------------------- B. screening for random term and background --------------------------
  M0 <- length(genotype0)

  for(i in 1:M0){
    attr(genotype0[[i]], "dimnames") <- NULL
  } # remove attr


  genotype0 <- lapply(1:M0, function(m)  replace(genotype0[[m]],is.na(genotype0[[m]]), 0  ) )
  # y for training without NA
  yc <- scale(y, center = T, scale = T)[-index_test]

  center <- attr(scale(y, center = T, scale = T),"scaled:center")
  scale <- attr(scale(y, center = T, scale = T),"scaled:scale")

  genotype <- genotype0
  #--------------------------- select the gene end------------------------------

  M <- length(genotype)
  regions_mat <- lapply(1:M, function(m) as.matrix(scale((genotype[[m]]), center= T, scale= F)) )

  tot_rare_common_var0 <- do.call("cbind", genotype)

  # center the big matrix with all snps
  tot_rare_common_var <- as.matrix(scale(tot_rare_common_var0, center= T, scale= F))

  #-------------------------- B. screening for fixed term --------------------------

  #
  tmp=lapply(1:M0, function(m){
    tmp1=genotype0[[m]]
    if(is.null(dim(tmp1))){
      tmp1=data.frame(tmp1)
    }
    tmp1
  })
  genotype0=tmp
  index_var_snp <- lapply(1:M0, function(m) which(getmafmat(genotype0[[m]]) > maf.filter.snp) )

  # obtain the common snps by using index

  var_snp <- lapply(1:M0, function(m) genotype0[[m]][, index_var_snp[[m]] ] )
  # var_snp <- lapply(1:M0, function(m) regions_mat[[m]][, index_var_snp[[m]] ] )

  # combine all snps from all genes
  tot_snp_only0 <- do.call("cbind", var_snp)

 # tot_snp_only0 <- readRDS(file = "/scale_wlg_persistent/filesets/project/nesi00464/test05062019/snp_condidate.rds")

  #X <- scale(genotype[-index_test,], center= T, scale=standardiseGenotypes_sd(genotype[-index_test,]))

  # screening ---------------------------------------
  Np_total <- dim(tot_snp_only0)[2]
  Np <- min(N, Np_total,10)
  # Best subset selection
  screening_fix <- bess.one(tot_snp_only0[-index_test,], yc, s = Np, family = "gaussian", normalize = TRUE)

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
   tot_snp_only <- as.matrix(scale(tot_snp_only0[,index_fixed_marker2], center= T, scale= T))
 # tot_snp_only <- tot_snp_only0

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
  # data = genotype[[2]]; N = N; type = weight.type;
  if(is.null(weight.type) == FALSE){

    if(M==1){
      w0 <- lapply(1:M, function(m) weight_function(genotype[[m]], N = N, type = weight.type ) )
    }else{
      w0 <- lapply(1:M, function(m) weight_function(genotype[[m]], N = N, type = "uw" ) )
      w0[1] <- lapply(1:1, function(m) weight_function(genotype[[m]], N = N, type = weight.type ) )
    }

    if(weight.type=="uw"){
      w  <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))) )
      pve <- 0.8
    }else if(weight.type=="beta"){
      w  <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))) )
      pve <- 0.8
    }else if(weight.type=="wss"){
      w  <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))) )
      pve <- 0.8
    }

    # data=genotype[[2]]; N = N; type = weight.type;

  }

  gc(verbose = FALSE)

#-------------------------- C. reconstruct the covarience matrix for random term--------------------------

  # 1.  Z original SNP matrix
  # ZList0 <- lapply(1:M, function(m)  dropNA(regions_mat[[m]])  )

  #  ZList0 <- lapply(1:M, function(m)  regions_mat[[m]] )

  NP <- lapply(1:M, function(m)  dim(regions_mat[[m]])[2] )
  ZList0 <- lapply(1:M, function(m)  scale(regions_mat[[m]], F, F) )

  # A_ZList_sd <- lapply(1:M, function(m) standardiseGenotypes_sd(ZList00[[m]]))  # get constant C which is variance of each SNPs

  # 2. matrix normalization

  # ZList0 <- lapply(1:M, function(m) scale(ZList00[[m]], center= T, scale=A_ZList_sd[[m]]) )

  #ZList0 <- ZList00

  # ZList0 <- ZList00 # for simulation data

  # 3. select the kernel function - Lin, RBF, IBS, POLY2
  ## weight function - uw, wss, beta
  ## calculate kernel matrix

  Cov_List <- list()

  for (m in 1:M0) {

    if(kernel[[m]] == "NA"){
      Cov_List[[m]] <- tcrossprod(ZList0[[m]])/NP[[m]]  # Linear effects for SNPs
    }else if(kernel[[m]]=="Lin" && is.null(weight.type) == TRUE){
      Cov_List[[m]] <- tcrossprod(ZList0[[m]])/NP[[m]]  # Linear effects for SNPs
    }else if(kernel[[m]]=="RBF" && is.null(weight.type) == FALSE){
      Cov_List[[m]]   <- rbfkernel(ZList0[[m]])/NP[[m]]  # RBF effects for SNPs
    }else if(kernel[[m]]=="POLY2" && is.null(weight.type) == FALSE){
      Cov_List[[m]] <- ((tcrossprod(ZList0[[m]]))**2)/NP[[m]]  # Poly 2 effects for SNPs
    }else if(kernel[[m]]=="Lin" && is.null(weight.type) == FALSE ){
      #  Cov_List[[m]] <- (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]]))/NP[[m]]
      Cov_List[[m]] <- (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]]))/NP[[m]]
      #  Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]])) )
    }else if(kernel[[m]]=="quadratic" && is.null(weight.type) == FALSE ){
      # Cov_List[[m]] <- ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2)/NP[[m]]
      Cov_List[[m]] <- ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2)/NP[[m]]
      #  Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2) )
    }else if(kernel=="NN" && is.null(weight.type) == FALSE ){
      Cov_List[[m]] <- NNnetwork(ZList0[[m]])/NP[[m]]
      #  Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2) )
    }else{
      stop ("currently only support Lin, RBF, POLY2 and IBS kernel and uw, beta, wss weight.type")
    }

  }

  # if(kernel=="Lin" && is.null(weight.type) == TRUE){
  #    Cov_List <- lapply(1:M, function(m) tcrossprod(ZList0[[m]]) )  # Linear effects for SNPs
  #  }else if(kernel=="RBF" && is.null(weight.type) == TRUE){
  #    Cov_List  <- lapply(1:M, function(m) (kernelMatrix(rbfdot(sigma = 0.3), ZList0[[m]])@.Data) ) # RBF effects for SNPs
  #  }else if(kernel=="POLY2" && is.null(weight.type) == TRUE){
  #    Cov_List <- lapply(1:M, function(m) (tcrossprod(ZList0[[m]]) + 1)**2 ) # Poly 2 effects for SNPs
  #  }else if(kernel=="Lin" && is.null(weight.type) == FALSE ){
  #    Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]]))/NP[[m]] )
  #  Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]])) )
  #  }else if(kernel=="quadratic" && is.null(weight.type) == FALSE ){
  #    Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2)/NP[[m]] )
  #  Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2) )
  #  }else{
  #    stop ("currently only support Lin, RBF, POLY2 and IBS kernel and uw, beta, wss weight.type")
  #  }

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

  # 5. add the pairwise interaction terms - random effects
  if(M0>1){ # if pairwise interaction exists

    pairs <- combinations(M0, 2, 1:M0)
    Npair <- dim(pairs)[1]
    Cov_int_pairs <- lapply(1:Npair, function(p) hadamard.prod(Cov_List[[pairs[p,1]]], Cov_List[[pairs[p,2]]]) ) # covariance of interaction pairs
    Cov_List[(M0+1):(M0+Npair)] <- Cov_int_pairs # add interaction kernel to main kernel list

    M <- M0 + Npair # add number of pairs to M

  }else{
    M <- M0
  }
 #  Cov_List[[m]] <- ((tcrossprod(ZList0[[m]]))**2)/NP[[m]]  NNnetwork(ZList0[[m]])/NP[[m]]
  if(M0==2){
    if(kernel[[1]]=="RBF"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
       Cov_List[[M+2]] <- ((tcrossprod(ZList0[[1]]))**2)/NP[[1]]
      Cov_List[[M+3]] <- tcrossprod(ZList0[[2]])/NP[[2]]
       Cov_List[[M+4]] <- ((tcrossprod(ZList0[[2]]))**2)/NP[[2]]
    }else if(kernel[[1]]=="NN"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]]   <-((tcrossprod(ZList0[[1]]))**2)/NP[[1]]
      Cov_List[[M+3]] <- tcrossprod(ZList0[[2]])/NP[[2]]
       Cov_List[[M+4]]   <- ((tcrossprod(ZList0[[2]]))**2)/NP[[2]]
    }else if(kernel[[1]]=="POLY2"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]]   <-NNnetwork(ZList0[[1]])/NP[[1]]
      Cov_List[[M+3]] <- tcrossprod(ZList0[[2]])/NP[[2]]
       Cov_List[[M+4]]   <- NNnetwork(ZList0[[2]])/NP[[2]]
    }
  }else if(M0==3){
    if(kernel[[1]]=="RBF"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]] <- ((tcrossprod(ZList0[[1]]))**2)/NP[[1]]
      Cov_List[[M+3]] <- tcrossprod(ZList0[[2]])/NP[[2]]
      Cov_List[[M+4]] <- ((tcrossprod(ZList0[[2]]))**2)/NP[[2]]
      Cov_List[[M+5]] <- tcrossprod(ZList0[[3]])/NP[[3]]
      Cov_List[[M+6]] <- ((tcrossprod(ZList0[[3]]))**2)/NP[[3]]
    }else if(kernel[[1]]=="NN"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]]   <- ((tcrossprod(ZList0[[1]]))**2)/NP[[1]]
      Cov_List[[M+3]] <- tcrossprod(ZList0[[2]])/NP[[2]]
      Cov_List[[M+4]]   <- ((tcrossprod(ZList0[[2]]))**2)/NP[[2]]
      Cov_List[[M+5]] <- tcrossprod(ZList0[[3]])/NP[[3]]
      Cov_List[[M+6]]   <- ((tcrossprod(ZList0[[3]]))**2)/NP[[3]]
    }else if(kernel[[1]]=="POLY2"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]]   <- NNnetwork(ZList0[[1]])/NP[[1]]
      Cov_List[[M+3]] <- tcrossprod(ZList0[[2]])/NP[[2]]
      Cov_List[[M+4]]   <- NNnetwork(ZList0[[2]])/NP[[2]]
      Cov_List[[M+5]] <- tcrossprod(ZList0[[3]])/NP[[3]]
      Cov_List[[M+6]]   <- NNnetwork(ZList0[[3]])/NP[[3]]
    }
  }else{
    if(kernel[[1]]=="RBF"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]] <- ((tcrossprod(ZList0[[1]]))**2)/NP[[1]]
    }else if(kernel[[1]]=="NN"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]]   <- ((tcrossprod(ZList0[[1]]))**2)/NP[[1]]
    }else if(kernel[[1]]=="POLY2"){
      Cov_List[[M+1]] <- tcrossprod(ZList0[[1]])/NP[[1]]
      Cov_List[[M+2]]   <- NNnetwork(ZList0[[1]])/NP[[1]]
    }
  }

  M <- length(Cov_List)

 # Cov_List <- lapply(1:M, function(m) StandardizeKernel(Cov_List[[m]]) )

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

  NP0 <- dim(tot_rare_common_var)[2]

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

  #  A      <- yc - rowSums(ZGammaU[,, drop=F]) - g

    # Beta

  #  s_beta <- solve( drop(((A_sigma2e0 + N/2)/B_sigma2e)) * hadamard.prod(crossprod(X), omega) +
  #                     drop((A_sigma2beta0 + N/2)/B_sigma2beta) * diag(1, Px)
  #  )

  #  m_beta <-  drop(((A_sigma2e0 + N/2)/B_sigma2e)) * s_beta %*% diag(E_w_beta) %*% t(X) %*% A

  # sigma2 for fixed effects

  #  A_sigma2beta <- (1/2) * N + A_sigma2beta0

  #  mtm_beta <-  crossprod(m_beta)

  #  tr_s_beta <- tr(s_beta)

  #  B_sigma2beta <- (1/2) * (mtm_beta + tr_s_beta) + B_sigma2beta0

    # B_sigma2beta <- B_sigma2beta0

    # Gamma w_beta

    # use X_E_w_beta[,-j] = X[,-j] %*% diag(E_w_beta[-j]) to replace (for computation purpose)
  #  X_E_w_beta <- X %*% diag(E_w_beta)

  #  logit_gamma_beta <- sapply(1:Px, function(j) logit(E_theta_beta[j]) -
  #                               (1/2) * drop((A_sigma2e0 + N/2)/B_sigma2e) * crossprod(X[,j]) %*% (crossprod(m_beta[j]) + s_beta[j,j]) +
  #                               drop(((A_sigma2e0 + N/2)/B_sigma2e)) * t(X[,j]) %*% ( A * m_beta[j] -
  #                                                                                       X_E_w_beta[,-j] %*% (m_beta[-j] * m_beta[j] + s_beta[-j,j])
  #                               )
  #  )

  #  E_w_beta <- sapply(1:Px, function(j) inv.logit(logit_gamma_beta[j])
  #  )

    # theta for fixed effects

  #  A_theta_beta <- sapply(1:Px, function(j) E_w_beta[j] + A_theta_beta0
  #  )

  #  B_theta_beta  <- sapply(1:Px, function(j) 1 - E_w_beta[j] + B_theta_beta0
  #  )

  #  E_theta_beta  <- sapply(1:Px, function(j) A_theta_beta[j]/(A_theta_beta[j] + B_theta_beta[j])
  #  )

    #  E_theta_beta <- rep(theta_beta, Px)

    # Omega
  #  omega  <- tcrossprod(E_w_beta) + hadamard.prod(diag(E_w_beta), (diag(1,Px)-diag(E_w_beta)))

  #  XGammaBeta <- X %*% diag(E_w_beta) %*% m_beta

    #----------------------# background #----------------------#

    B0      <- yc - XGammaBeta - rowSums(ZGammaU[,, drop=F])
  #  B0      <- yc  - rowSums(ZGammaU[,, drop=F])


    # Beta

    s_b <- solve( drop(((A_sigma2e0 + N/2)/B_sigma2e)) * tglg +
                    drop((A_sigma2b0 + N/2)/B_sigma2b) * diag(1, DIM0)
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
  #  CtC4  <- tr( hadamard.prod(crossprod(X), omega) %*% s_beta
  #  )

    CtC5  <- sapply(1:M, function(m) tr( drop(E_w_u[[m]]) * crossprod(Z[[m]]) %*% s_u[[m]]
    )
    )

    CtC6  <- tr(tglg %*% s_b
    )

    B_sigma2e <- B_sigma2e0 + (1/2) * ( crossprod(yc - CtC1 - rowSums(CtC2) - CtC3 ) +
                                      #    sum(CtC4) +
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

   # lb_pbeta <- N/2 * log(2 * pi) - 1/2 * (A_sigma2beta/B_sigma2beta) * (crossprod(m_beta) + tr(s_beta))

  #  lb_psigma2_beta <- A_sigma2beta0 * log(B_sigma2beta0) - lgamma(A_sigma2beta0) -
  #    (A_sigma2beta0 + 1) * (B_sigma2beta/(A_sigma2beta - 1)) + B_sigma2beta0 * (A_sigma2beta/B_sigma2beta)

   # lb_pgamma_beta <- sapply(1:Px, function(j) E_w_beta[j] * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
   #                            (1 - E_w_beta[j]) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
   # )

  #  lb_ptheta_beta <- sapply(1:Px, function(j) (A_theta_beta0 - 1) * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
   #                            (B_theta_beta0 - 1) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
   # )

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

  #  lb_qsigma2_beta <- - A_sigma2beta - log(B_sigma2beta) * lgamma(A_sigma2beta) +
  #    (1 + A_sigma2beta) * digamma(A_sigma2beta)

    #  lb_qbeta <- - 1/2 * log(exp(determinant(s_beta)$modulus[1]) + 1e-10) - 1/2 * log(2 * pi)

   # lb_qbeta <- - 1/2 * determinant(s_beta)$modulus[1] - 1/2 * log(2 * pi)

  #  lb_qgamma_beta <- sapply(1:Px, function(j) E_w_beta[j] * log(E_w_beta[j] + 1e-10) +
  #                             (1 - E_w_beta[j]) * log(1 - E_w_beta[j] + 1e-10)
  #  )

   # lb_qtheta_beta <- sapply(1:Px, function(j) lgamma(A_theta_beta[j] + B_theta_beta[j]) - lgamma(A_theta_beta[j]) - lgamma(B_theta_beta[j]) +
   #                            (A_theta_beta[j] - 1) * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
   #                            (B_theta_beta[j] - 1) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
   # )

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
     # lb_psigma2_beta + sum(lb_pgamma_beta) + sum(lb_ptheta_beta) +
      sum(lb_pu) + sum(lb_psigma2_u) + sum(lb_pgamma_u) + sum(lb_ptheta_u) +
      lb_pb + lb_psigma2_b +
      lb_psigma2_e  -
    #  lb_qsigma2_beta - lb_qbeta - sum(lb_qgamma_beta) - sum(lb_qtheta_beta) -
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
    if (i == max_iter ) {message("[2] VB did not converge!\n")
      }

  }

  if (log_delta < epsilon_conv) {
    obj <- structure(list(X = X, Z = Z, genotype = tot_snp_only0, genotype0 = tot_snp_only, gamma_lambda = gamma_lambda,
                          B_sigma2beta= B_sigma2beta, index_test = index_test, B_sigma2e = B_sigma2e, B_sigma2u = B_sigma2u,
                          post_sigma2b =post_sigma2b, post_sigma2u = post_sigma2u, post_sigma2beta = post_sigma2beta,
                          post_sigma2e = post_sigma2e,N = N, M = M, Px = Px, DIM = DIM, DIM0 = DIM0, Delta = log_delta,
                          E_theta_u = E_theta_u, E_theta_beta = E_theta_beta, m_b = m_b, s_b = s_b,
                          Cov_List = Cov_List, m_beta = m_beta,r = M,m_u = m_u,s_u = s_u, y=yc, scale = scale, center = center,
                          E_w_beta = E_w_beta, E_w_u = E_w_u, Px =Px,index_fixed_marker = index_fixed_marker2,
                          L = L[2:i]), UD = UListD, class = "vb")
    return(obj)

  }

}

