# STEP A1 (generate the input data matrix)
# Modified #

#' @import MASS
STEP_A1 <- function(gene_path, exp_path, meth_path,readFun=readRDS){
  if(length(gene_path)!=length(exp_path)){
    stop("The number of genes listed in the genotype data is not the same as that in gene expression. Note for missing, list NA.")
  }
  if(length(gene_path)!=length(meth_path)){
    stop("The number of genes listed in the genotype data is not the same as that in methylation. Note for missing, list NA.")
  }
  if(length(exp_path)!=length(meth_path)){
    stop("The number of genes listed in the methylation data is not the same as that in gene expression. Note for missing, list NA.")
  }

  n1 <- length(gene_path)
  gene_file <- gene_path
  exp_file <- exp_path
  meth_file <- meth_path

  comb_all <- list()

  for (i in 1:n1) {
    comb_all[[i]]=list()
    if(!is.na(gene_file[i])){
      gene_RDS <- readFun(gene_file[i])
      comb_all[[i]]=append(comb_all[[i]],list(gene_RDS))
    }
    if(!is.na(exp_file[i])){
      exp_RDS <- readFun(exp_file[i])
      comb_all[[i]]=append(comb_all[[i]],list(exp_RDS))
    }
    if(!is.na(meth_path[i])){
      meth_RDS <- readFun(meth_path[i])
      comb_all[[i]]=append(comb_all[[i]],list(meth_RDS))
    }
  }

  gene_file=gene_file[!is.na(gene_file)]
  if(length(gene_file)==0) stop('Current implementation expect to have genotype files!')
  geno_all <- lapply(1:n1, function(i) readFun(gene_file[i]))
  return(list(comb_all = comb_all,geno_all = geno_all))
}


