# Read file - input can only be tab-delimited text or csv
read_file <- function(filename){
  file_type <- gsub(".*\\.","", filename)
  if(file_type == "txt"){
    read.delim(filename, row.names = 1)
  }
  if(file_type == "txt" | file_type == "tsv"){
    read.csv(filename, row.names = 1)
  }
}


# param - z_score_row, log2
transform_data <- function(df, param){
  # Compute z-scores
  if(is.null(param))
    return(df) 
  
  # row
  if("z_row" %in% param)
    df <- t(z_score(df, 1)) # apply to row and transpose back to original dimensions
    
  # column
  if("z_col" %in% param)
    df <- z_score(df, 2)
  
  # Compute log2
  if("log2" %in% param)
    df <- log2(df + 1)
  
  return(df) 
}


# Z score formula
# row - apply by row (default)

z_score <- function(df, index=1){
  apply(df, index, function(x) (x - mean(x)) / sd(x))
}

# Return FALSE if NULL
check <- function(param, value){
  if(is.null(param)) return(F)
  
  ifelse(value %in% param, T, F)
}