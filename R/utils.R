#' Round a decimal number to two digits
#'
#' @param x a number.
#'
#' @return A string of the input number rounded to two digits.
rd = function(x) {
  return(sprintf("%.2f", x))
}


#' Round a p-value to three digits
#'
#' @param x a number.
#'
#' @return A string of the input number rounded to three digits or string
#' containing "<0.001" if the number is less than 0.001.
rd_p = function(x){
  if(x<0.001){
    return("<0.001")
  }
  else{
    return(paste0(" ", sprintf("%.3f", x)))
  }
}

#' Round a percentage to two digits
#'
#' @param x a number.
#'
#' @return A string of the input number as percentage rounded to two digits
#' (i.e. 0 digits in terms of percentage).
rd_pct = function(x) {
  return(sprintf("%.0f", x*100))
}
