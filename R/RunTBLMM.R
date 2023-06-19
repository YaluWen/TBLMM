#' Run TBLMM and make prediction
#' @importFrom utils read.table
#' @param phenofile path to the phenotype files that have 3 columns (PID, IID and Y). The individuals to be predicted should be marked as \code{NA}.
#' @param gene_path path to each genotype file (current implementation expects to have genotype files).
#' @param exp_path path to each gene expression file.
#' @param meth_path path to each methylation file.
#' @param kernel kernel function for each omics (default is linear kernel.)
#' @param readFun the function used to read gene expression, genotype and methylation files. The default is \code{readRDS}.
#' @param verbose print the fitting process (default is \code{TRUE}).
#' @param modelFit return all model fits including parameter estimates (default is \code{FALSE}, and only predicted values to be exported).
#'
#' @return A list of predicted values.
#' @export
#'
#' @examples
#' datadir= system.file("extdata", package = "TBLMM")
#' gene_path=paste0(datadir,"/geno",1:2,".RDS")
#' expr_path=paste0(datadir,"/expr",1:2,".RDS")
#' meth_path=paste0(datadir,"/meth",1:2,".RDS")
#' phenofile=paste0(datadir,"/pheno.txt")
#' test=RunTBLMM(phenofile=phenofile,gene_path=gene_path,exp_path=expr_path, meth_path=meth_path,
#' kernel=c("linear","linear","linear"),readFun=readRDS)
RunTBLMM<-function(phenofile,gene_path=NULL,exp_path=NULL, meth_path=NULL,kernel=NULL,readFun=readRDS,
                   verbose=TRUE, modelFit=FALSE){
  maf_f = 0.001
  print(getwd())
  if(is.null(gene_path) & is.null(exp_path) & is.null(meth_path)){
    stop("Must include some data!")
  }
  if(is.null(kernel)){
    cat("No kernels have been specified, and only linear kernel will be used \n")
    kernel=c("linear","linear","linear")
  }

  # 1 STEP A1 - get data matrix
  cat("[STEP 0-A1] Data matrix:","\n")
  mat <- STEP_A1(gene_path, exp_path, meth_path)
  nm <- length(mat$comb_all)
  cat("[STEP 0-A1] Done!","\n")

  # 1 STEP B1 - generate the kernel list for single gene analysis
  kernel[kernel=="linear"]="Lin"
  kernel[kernel=="poly2"]="POLY2"
  kernel[kernel=="neural"]="NN"
  kernel=paste0(kernel,collapse = '.')
  cat("[STEP 0-B1] Lernel list:","\n")
  kernel_list <- STEP_B1(mat$comb_all, kernel)
  cat("[STEP 0-B1] Done!","\n")

  # 2 Get phenotypes #
  cat("[STEP 0-C1] Phenotypes:","\n")
  y_na <- read.table(phenofile, header = T)
  if(ncol(y_na)!=3) stop("Phenotype must have 3 columns with column name being PID,IID and Y")
  if(sum(colnames(y_na)!=c("PID","IID","Y"))!=0) stop("Phenotype must have 3 columns with column name being PID,IID and Y")
  if(sum(is.na(y_na))==0) stop("No individuals with missing phenotypes, and no prediciton will be performed.")
  if(nrow(y_na)!=nrow(mat$comb_all[[1]][[1]])) stop("The number of individuals in omics is different from phenotypes")
  IDpred=y_na[is.na(y_na$Y),1:2]
  y_na=y_na$Y
  cat("[STEP 0-C1] Done!","\n")


  # 3 Single gene analyais
  cat("[STEP 1] Single gene analysis:","\n")
  fit1 <- list()
  for (s1 in 1:nm) {
    cat("[STEP 1] Gene:",s1,"\n")
    fit1[[s1]] <- vb_fit_STEP1(y = y_na, genotype0 = mat$comb_all[[s1]], max_iter = 25000, weight.type = "uw",
                               maf.filter.snp = maf_f, epsilon_conv = 1e-4, verbose = verbose, kernel = kernel_list[[s1]])
  }
  cat("[STEP 1] Done!","\n")

  # 4 STEP C2 - gen the composite kernels for the STEP 2
  cat("[STEP C2] Composite kernels:","\n")
  cov_kernel <- STEP_C2(fit1, ng = length(fit1))
  cat("[STEP C2] Done!","\n")

  cat("[STEP 2] BLMM:","\n")
  fit2 <- vb_fit_STEP2(y = y_na, mat_add = cov_kernel$mat_add, kernel_ran = cov_kernel$cmat_kernel, genoback = mat$geno_all,
                       maf.filter.snp = 0.001, max_iter = 25000, weight.type = "uw",  epsilon_conv = 1e-4, Bsigma2beta = 1,
                       theta_beta = 0.1, theta_u = 0.1, verbose = TRUE)
  cat("[STEP 2] Done!","\n")

  pred =vb_predictive(fit2, maf_beta = 0.5, maf_u = 0.5, B = "Yes")
  pred$ID=IDpred
  if(!modelFit){
    return(cbind(IDpred,pred$y_hat))
  }
  if(modelFit){
    return(pred)
  }
}
