#'
#' @name print.bcgfrailev
#' @title Print bcgfrailev model
#' @description Generics to print the S3 class bcgfrailev.
#' @details Calls \code{print.bcfrailph()}.
#' 
#' @param x A class \code{bcgfrailev} object. 
#' @param ... ignored
#'  
#' @export 
#  (deprecated) @S3method print bcgfrailev
#' 

print.bcgfrailev<- function( x, ... ) {
  
  if (!inherits(x,"bcgfrailev")) stop("Class of the argument must be bcgfrailev.")

  NextMethod() # calls print.bcfrailph()
}

