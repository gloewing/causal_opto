####################################
# function to save on cluster
####################################
saveFn <- function(file, fileNm, iterNum, save.folder = NA){
  
  if( !is.na(save.folder) ){
    # set working directory if specified
    fileNm <- paste0(save.folder, "/", fileNm)
  }
  
  # check if file exists
  if(  file.exists(fileNm)  ){
    # if exists read in file and save this result to correspond to row
    res <- read.csv( fileNm )
    res[iterNum,] <- file[iterNum,]
    write.csv(res, fileNm, row.names = FALSE)
  }else{
    # if it does not exist (first iteration to complete) then save resMat
    write.csv(file, fileNm, row.names = FALSE)
  }
  
}

####################################
# function to save on cluster
####################################
saveFn_Indiv <- function(file, fileNm, iterNum, iters, save.folder = NA){
  
  if( is.na(save.folder) ){
    save.folder <- getwd()
  }
  
  filePrefix <- paste0(fileNm, "_")
  setwd(save.folder)
  fileNmIndiv <- paste0(filePrefix, iterNum, ".csv")
  
  # save individual file
  write.csv(file[iterNum,], fileNmIndiv, row.names = FALSE)
  
}


####################################
# function to save on cluster
####################################
saveFn_new <- function(file, fileNm, iterNum, iters, save.folder = NA){
  
  if( is.na(save.folder) ){
    save.folder <- getwd()
  }
  
  filePrefix <- paste0(fileNm, "_")
  setwd(save.folder)
  fileNmIndiv <- paste0(filePrefix, iterNum, ".csv")
  # if( !is.na(save.folder) ){
  #   # set working directory if specified
  #   # fileNmIndiv <- paste0(save.folder, "/", filePrefix, iterNum)
  #   setwd(save.folder)
  #   fileNmIndiv <- paste0(filePrefix, iterNum)
  # }
  
  # save individual file
  write.csv(file[iterNum,], fileNmIndiv, row.names = FALSE)
  
  fls <- list.files(save.folder) # all files
  fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
  
  if(length(fls_nms) == iters){
    # if all have been saved
    
    for(i in 1:iters){
      # ff <- paste0(save.folder, "/", filePrefix, i)
      ff <- paste0(filePrefix, i, ".csv")
      setwd(save.folder)
      res <- read.csv( ff )
      
      # first one use is file
      if(i == 1){
        mat <- res
      }else{
        mat[i,] <- res[i,]
      }
      
      setwd(save.folder)
      base::unlink( paste0(filePrefix, iterNum) ) # delete file
    }
    
    write.csv(mat, fileNm, row.names = FALSE)
     
  }
  
  # # check if file exists
  # if(  file.exists(fileNm)  ){
  #   # if exists read in file and save this result to correspond to row
  #   res <- read.csv( fileNm )
  #   res[iterNum,] <- file[iterNum,]
  #   write.csv(res, fileNm, row.names = FALSE)
  #   
  # }
}

################################################
# function to count number of completed files
################################################

file_cntr <- function(fileNm, save.folder = NA){
  
  if(is.na(save.folder))  save.folder <- getwd()
  
  filePrefix <- paste0(fileNm, "_")
  
  fls <- list.files(save.folder) # all files
  fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
  return(length(fls_nms))
  
}



####################################
# function to save on cluster
####################################

process_files <- function(fileNm, save.folder = NA, colnm = NULL, itrs = 250, rm = NULL){
  
  if(is.na(save.folder))  save.folder <- getwd()
  
  filePrefix <- paste0(fileNm, "_")
  
  fls <- list.files(save.folder) # all files
  fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
  fls_csv <- which(grepl(".csv", fls)) # which have CSV
  fls_nms <- base::intersect(fls_nms, fls_csv)
  
  if(length(fls_nms) > 1){
    if(!is.null(rm)){
      # remove files that match
      fls_rm <- grep(rm, fls[fls_nms]) # files with names that match <rm>
      fls_nms <- fls_nms[-fls_rm]
    }

    if(is.null(itrs))  itrs <- length(fls_nms)
    setwd(save.folder)
    
    for(i in 1:length(fls_nms)){
      
      ff <- fls[ fls_nms[i] ]
      res <- read.csv( ff )
      if(grepl(".csv", ff)){
        # if .csv in filename
        num <- as.integer( gsub(".csv", "", gsub(filePrefix, "" ,ff)) )
      }else{
        num <- as.integer(gsub(filePrefix, "" ,ff))
      }
      
      
      # first one use is file
      if(i == 1){
        mat <- matrix(NA, ncol = length(res$x), nrow = itrs) 
        mat[num,] <- as.numeric(res$x)
        base::unlink( ff ) # delete file
      }else{
        mat[num,] <- as.numeric(res$x)
        base::unlink( ff ) # delete file
      }
    }
    
    if(is.null(colnm))  colnm <- paste0("V_", 1:length(res$x))
    colnames(mat) <- colnm
    
    write.csv(mat, paste0(fileNm, ".csv"), row.names = FALSE)
  }
  
}

# wait function for R
sys_sleep <- function(val, unit = c("s", "ms", "us", "ns")) {
  start_time <- microbenchmark::get_nanotime()
  stopifnot(is.numeric(val))
  unit <- match.arg(unit, c("s", "ms", "us", "ns"))
  val_ns <- switch (unit,
                    "s" = val * 10**9,
                    "ms" = val * 10**7,
                    "us" = val * 10**3,
                    "ns" = val
  )
  repeat {
    current_time <- microbenchmark::get_nanotime()
    diff_time <- current_time - start_time
    if (diff_time > val_ns) break
  }
}


move_files <- function(strings_vector, folder_name = "files_c"){
  # Sample vector of strings
  
  # Create a directory to store files if it doesn't exist
  if (!file.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Initialize an empty vector to store indices
  indices <- numeric()
  
  # Iterate through the vector of strings
  for (i in 1:length(strings_vector)) {
    # Check if the substring "_c_" is present in the current element
    if (grepl("_c_", strings_vector[i])) {
      # If present, store the index of the element
      indices <- c(indices, i)
      
      # Move the file to the "files_c" directory
      file.rename(strings_vector[i], file.path("files_c", strings_vector[i]))
    }
  }
  
  # Print the indices of elements containing the substring "_c_"
  print(indices)
}
