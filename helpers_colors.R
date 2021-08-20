require(RColorBrewer)
# Example of function use:
# library(dplyr)
# get_col_pal("RdBu") %>%
#   get_col_grad(n=3) %>%
#   scales::show_col()


# Get color palette by specifying: 
#' brew_pal - character or NA. exact name of RColorBrewer palette
#'OR custom - NA or character vector with colors that will create a gradient
# rev = logical. reverse color palette? 
get_col_pal<- function(brew_pal, custom=NA, rev=F){
  # Get custom palette based on vector of colors
  if(!is.na(custom)) {
    colorRampPalette(custom)    #c(low_col, high_col))
  }
  
  if(!is.na(brew_pal)){
    # Get max colors of palette
    max_n <- brewer.pal.info$maxcolors[grep(brew_pal, rownames(brewer.pal.info))]
    # Create brewer pal
    # Reverse if required
    if(rev){
      brew <- rev(brewer.pal(n = max_n, name = brew_pal))
    }else{
      brew <- brewer.pal(n = max_n, name = brew_pal)
    }
    # Return palette
    colorRampPalette(brew)
  }
}

# Get color gradient by specifying:
# pal = character vector, names of colors (hex)
# n = numeric. number of colors in gradient
get_col_grad <- function(pal, n){
  # Create and return gradient
  pal(n)
}


# Get colours for a vector
#' v = character vector. unique elements
#' colRamp = colorRampPalette with colors of palette
#' rearr = logical. rearrange vector so it's in random order?
get_element_colors <- function(v, colRamp, rearr = F){
  require(RColorBrewer)
  # Get unique elements and rearrange
  v <- unique(v)
  
  # Rearrange   
  if(rearr)
    v <- sample(v)
  
  # Get color gradient
  myColors <- colRamp(length(v))
  names(myColors) <- v
  return(myColors)
}