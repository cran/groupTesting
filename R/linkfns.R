#' Link Functions in the Class of Generalized Linear Models
#'
#' This function provides characteristics of common link functions (logit, probit, and complementary log-log). Specifically, based on the link name, the inverse link function, its first and second derivatives, the original link function, and the first derivative of the original link are returned.
#'
#' @param fn.name One of the three: "logit", "probit", and "cloglog".
#'
#' @return A list with components:
#' \item{g}{The inverse link function corresponding to "logit", "probit", or "cloglog".}
#' \item{dg}{The first derivative of g.}
#' \item{d2g}{The second derivative of g.}
#' \item{gInv}{The original link function.}
#' \item{dgInv}{The first derivative of gInv.}
#' 
#' @export
#' 
#' @examples
#'
#' ## Logit link functions
#' # Try: 
#' glmLink("logit")
#'
#' link <- glmLink("logit")
#' ## Evaluate inverse link and its first derivative
#' link$g(0)
#' link$dg(0)
#'
#' ## Evaluate the original link
#' link$gInv(0.5)
#' 
glmLink <- function(fn.name=c("logit", "probit", "cloglog")){
  fn.name <- match.arg(fn.name)
  if(fn.name == "logit"){
    g <- function(u){
      stats::plogis(u)
    }
    dg <- function(u){
      p <- stats::plogis(u)
      p*(1-p)
    }
    d2g <- function(u){
      p <- stats::plogis(u)
      p*(1-p)*(1-2*p)
    }
    g.inv <- function(u){
      stats::qlogis(u)
    }
    dg.inv <- function(u){
      1/(u*(1-u))
    }
  }
  if(fn.name == "probit"){
    g <- function(u){
      stats::pnorm(u)
    }
    dg <- function(u){
      exp(-u^2/2)/sqrt(2*pi)
    }
    d2g <- function(u){
      -u*exp(-u^2/2)/sqrt(2*pi)
    }
    g.inv <- function(u){
      stats::qnorm(u)
    }
    dg.inv <- function(u){
      1/stats::dnorm(stats::qnorm(u))
    }
  }
  if(fn.name == "cloglog"){
    g <- function(u){
      1 - exp(-exp(u))
    }
    dg <- function(u){
      e.u <- exp(u)
      e.u*exp(-e.u)
    }
    d2g <- function(u){
      e.u <- exp(u)
      e.u*exp(-e.u)*(1-e.u)
    }
    g.inv <- function(u){
      log(-log(1-u))
    }
    dg.inv <- function(u){
      -1/((1-u)*log(1-u))
    }
  }
  list(
    "g" = g,
    "dg" = dg,
    "d2g" = d2g,
    "gInv" = g.inv,
    "dgInv" = dg.inv
  )
}
