# STEP B1 (generate kernel list)

STEP_B1 <- function(mat, kernel){
  
  nm <- length(mat)
  int_type_list <- strsplit(kernel,'.',fixed=TRUE)
  
  kernel_list <- list()
  
  for (i in 1:nm) {
    np <- length(mat[[i]])
    if(np==1){
      kernel_list[[i]] <- int_type_list[[1]][1]
    }else if(np==2){
      kernel_list[[i]] <- int_type_list[[1]][1:2]      
    }else if(np==3){
      kernel_list[[i]] <- int_type_list[[1]][1:3]      
    }
  }

  return(kernel_list)
}


