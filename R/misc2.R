# ----------------------------------------------------------------------
# Additional auxiliary functions. 
# Could be merged into misc.R after testing
# ----------------------------------------------------------------------

# Return true 
is.scalar <- function(x) {
  if (!(length(x) == 1 &
        is.numeric(x)) &
        !(is.na(x))) {
    stop("Input argument is not a scalar")
  } else {
    return(TRUE)
  }
}
