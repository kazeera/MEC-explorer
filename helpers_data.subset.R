# Subset table by cols_of_interest that match the list_to_match
get_to_keep <- function(all_elements, list_to_match){
  all_elements %in% list_to_match
}

# Subset to the indices specificed by rows and cols arguments 
subset_table <- function(df, rows = NULL, cols = NULL){
  # keep all rows/cols if not indicated
  if(is.null(rows))
    rows <- 1:nrow(df)

  if(is.null(cols))
    cols <- 1:ncol(df)
  
  # subset and return
  df[rows,cols,drop=F]
}