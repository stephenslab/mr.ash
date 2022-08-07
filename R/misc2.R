# ----------------------------------------------------------------------
# Additional auxiliary functions.
# Could be merged into misc.R after testing
# ----------------------------------------------------------------------

# Return true
is.scalar <- function(x) {
  if (!(length(x) == 1 &&
        is.numeric(x)) &&
        !(is.na(x))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.int <- function(x) {

  x %% 1 == 0

}
