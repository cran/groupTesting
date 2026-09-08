#' EM Algorithm to Estimate the Prevalence of a Disease from Group Testing Data
#'
#' This function implements an expectation-maximization (EM) algorithm to find the maximum likelihood estimate (MLE) of a disease prevalence, p, based on group testing data. The EM algorithm, which is outlined in Warasi (2023), can model pooled testing data arising from a wide range of group testing protocols used in practice, including hierarchical and array testing (Kim et al., 2007).
#'
#' @useDynLib groupTesting, .registration=TRUE
#' 
#' @param p0 Initial value of the prevalence \eqn{p}, between 0 and 1.
#' @param gtData A matrix or data.frame consisting of the pooled test outcomes and other information from a group testing application. Must be structured according to the description and example provided in \code{\link{gtData}}.
#' @param covariance Logical. If TRUE, the variance is calculated at the MLE.
#' @param nburn The number of initial Gibbs iterates to be discarded.
#' @param ngit The number of Gibbs iterates to be used in the E-step after discarding the initial \code{nburn} iterates as burn-in.
#' @param maxit The maximum number of EM steps (iterations) allowed in the EM algorithm.
#' @param tol Convergence tolerance used in the EM algorithm.
#' @param tracing When TRUE, progress in the EM algorithm is displayed.
#' @param conf.level Confidence level for the Wald confidence interval.
#'
#' @details
#'
#' The EM algorithm consists of an E-step and an M-step that are 
#' performed alternately. In the E-step, the expectation of the
#' complete-data log-likelihood is approximated, whereas in the
#' M-step, this quantity is maximized with respect to the prevalence 
#' \eqn{p} using a closed-form expression. The EM iterations continue
#' until convergence of the prevalence estimates is
#' achieved; see Warasi (2023) for further details.
#' 
#' The variance of the MLE is calculated using the missing data principle and the method outlined in Louis (1982).
#' 
#' @return A list with components:
#' \item{param}{The MLE of the disease prevalence.}
#' \item{covariance}{Estimated variance for the disease prevalence.}
#' \item{iterUsed}{The number of EM iterations performed.}
#' \item{convergence}{0 if the EM algorithm converges successfully and 1 if the iteration limit \code{maxit} has been reached.}
#' \item{summary}{Estimation summary with Wald confidence interval.}
#' 
#' @export
#'
#' @references
#' Kim HY, Hudgens M, Dreyfuss J, Westreich D, and Pilcher C. (2007). Comparison of Group Testing Algorithms for Case Identification in the Presence of Testing Error. \emph{Biometrics}, 63:1152-1163.
#' 
#' Litvak E, Tu X, and Pagano M. (1994). Screening for the Presence of a Disease by Pooling Sera Samples. \emph{Journal of the American Statistical Association}, 89:424-434.
#' 
#' Liu A, Liu C, Zhang Z, and Albert P. (2012). Optimality of Group Testing in the Presence of Misclassification. \emph{Biometrika}, 99:245-251.
#' 
#' Louis T. (1982). Finding the Observed Information Matrix when Using the EM algorithm. \emph{Journal of the Royal Statistical Society: Series B}, 44:226-233.
#' 
#' Warasi M. (2023). groupTesting: An R Package for Group Testing Estimation. \emph{Communications in Statistics-Simulation and Computation}, 52:6210-6224.
#' 
#' @seealso
#' \code{\link{hier.gt.simulation}} and \code{\link{array.gt.simulation}} for group testing data simulation, \code{\link{glm.gt}} for group testing regression models, and \code{\link{gtData}} for information about the required data structure.
#'
#' @examples
#' 
#' library(groupTesting)
#' 
#' ## To illustrate 'prop.gt', we use data simulated by 
#' ## the R functions 'hier.gt.simulation' and 'array.gt.simulation'.
#' 
#' ## The simulated data-structures are consistent  
#' ## with the data-structure required for 'gtData'.
#' 
#' ## Example 1: MLE from 3-stage hierarchical group testing data.
#' ## The data used is simulated by 'hier.gt.simulation'.
#'
#' N <- 90               # Sample size
#' S <- 3                # 3-stage hierarchical testing
#' psz <- c(6,2,1)       # Pool sizes used in stages 1-3
#' Se <- c(.95,.95,.98)  # Sensitivities in stages 1-3
#' Sp <- c(.95,.98,.96)  # Specificities in stages 1-3
#' assayID <- c(1,2,3)   # Assays used in stages 1-3
#' p.t <- 0.05           # The TRUE parameter to be estimated
#'
#' # Simulating data:
#' set.seed(123)
#' gtOut <- hier.gt.simulation(N,p.t,S,psz,Se,Sp,assayID)$gtData
#'
#' # EM algorithm:
#' pStart <- p.t + 0.2   # Initial value
#' res <- prop.gt(p0=pStart,gtData=gtOut,covariance=TRUE,
#'                nburn=2000,ngit=5000,maxit=200,tol=1e-03,
#'                tracing=TRUE,conf.level=0.95)
#' 
#' # Estimation results:
#' # >  res
#' # $param
#' # [1] 0.05158
#' 
#' # $covariance
#' #              [,1]
#' # [1,] 0.0006374296
#' 
#' # $iterUsed
#' # [1] 4
#' 
#' # $convergence
#' # [1] 0
#' 
#' # $summary
#' #      Estimate Std.Err 95%lower 95%upper
#' # prop    0.052   0.025    0.002    0.101
#' 
#' ## Example 2: MLE from two-dimensional array testing data.
#' ## The data used is simulated by 'array.gt.simulation'.
#' 
#' N <- 100             # Sample size
#' protocol <- "A2"     # 2-stage array without testing the initial master pool
#' n <- 5               # Row/column size
#' Se <- c(0.95, 0.95)  # Sensitivities
#' Sp <- c(0.98, 0.98)  # Specificities
#' assayID <- c(1, 1)   # The same assay in both stages
#' p.true <- 0.05       # The TRUE parameter to be estimated
#' 
#' # Simulating data:
#' set.seed(123)
#' gtOut <- array.gt.simulation(N,p.true,protocol,n,Se,Sp,assayID)$gtData
#' 
#' # Fit the model:
#' pStart <- p.true + 0.2  # Initial value
#' res <- prop.gt(p0=pStart,gtData=gtOut,covariance=TRUE)
#' print(res)
#'
#' \donttest{
#' ## Example 3: MLE from non-overlapping initial pooled responses.
#' ## The data used is simulated by 'hier.gt.simulation'. 
#' 
#' ## Note: With initial pooled responses, our MLE is equivalent  
#' ## to the MLE in Litvak et al. (1994) and Liu et al. (2012).
#' 
#' N <- 1000           
#' psz <- 5             
#' S <- 1             
#' Se <- 0.95           
#' Sp <- 0.99           
#' assayID <- 1       
#' p.true <- 0.05      
#' 
#' set.seed(123)
#' gtOut <- hier.gt.simulation(N,p.true,S,psz,Se,Sp,assayID)$gtData
#' 
#' pStart <- p.true + 0.2  
#' res <- prop.gt(p0=pStart,gtData=gtOut,
#'                covariance=TRUE,nburn=2000,ngit=5000,
#'                maxit=200,tol=1e-03,tracing=TRUE)
#' print(res)
#'
#' ## Example 4: MLE from individual (one-by-one) testing data.
#' ## The data used is simulated by 'hier.gt.simulation'.
#' 
#' N <- 1000            
#' psz <- 1             
#' S <- 1               
#' Se <- 0.95         
#' Sp <- 0.99          
#' assayID <- 1         
#' p.true <- 0.05      
#' 
#' set.seed(123)
#' gtOut <- hier.gt.simulation(N,p.true,S,psz,Se,Sp,assayID)$gtData
#'
#' pStart <- p.true + 0.2   
#' res <- prop.gt(p0=pStart,gtData=gtOut,
#'                covariance=TRUE,nburn=2000,
#'                ngit=5000,maxit=200,
#'                tol=1e-03,tracing=TRUE)
#' print(res)
#' 
#' ## Example 5: Using pooled testing data.
#'
#' # Pooled test outcomes:
#' Z <- c(1, 0, 1, 0, 1, 0, 1, 0, 0)     
#' 
#' # Pool sizes used:
#' psz <- c(6, 6, 2, 2, 2, 1, 1, 1, 1)
#' 
#' # Pool-specific Se & Sp:
#' Se <- c(.90, .90, .95, .95, .95, .92, .92, .92, .92)
#' Sp <- c(.92, .92, .96, .96, .96, .90, .90, .90, .90)
#' 
#' # Assays used:
#' Assay <- c(1, 1, 2, 2, 2, 3, 3, 3, 3)
#' 
#' # Pool members:
#' Memb <- rbind( 
#'    c(1, 2,  3,  4,  5,  6),
#'    c(7, 8,  9, 10, 11, 12),
#'    c(1, 2, -9, -9, -9, -9),
#'    c(3, 4, -9, -9, -9, -9),
#'    c(5, 6, -9, -9, -9, -9),
#'    c(1,-9, -9, -9, -9, -9),
#'    c(2,-9, -9, -9, -9, -9),
#'    c(5,-9, -9, -9, -9, -9),
#'    c(6,-9, -9, -9, -9, -9)
#' )
#' # The data-structure suited for 'gtData':
#' gtOut <- cbind(Z, psz, Se, Sp, Assay, Memb)
#' 
#' # Fitting the model:
#' pStart <- 0.10
#' res <- prop.gt(p0=pStart,gtData=gtOut,
#'                covariance=TRUE,nburn=2000,
#'                ngit=5000,maxit=200,
#'                tol=1e-03,tracing=TRUE)
#' print(res)
#' }
#'
prop.gt <- function(p0,gtData,covariance=FALSE,nburn=2000,ngit=5000,maxit=200,tol=1e-03,tracing=TRUE,conf.level=0.95){

  # Validating the initial prevalence value p0
  if(!is.numeric(p0) || length(p0) != 1L || !is.finite(p0)) stop("p0 must be a finite numeric value.")
  if(p0 < 0 || p0 > 1) stop("p0 must be between 0 and 1.")
  p0 <- pmin(pmax(p0, 1e-12), 1-1e-12)

  # Validating the group testing data object gtData
  gtData <- as.matrix(gtData)
  if(!is.numeric(gtData)) stop("gtData must have numeric values.")
  if(ncol(gtData) < 6) stop("gtData must have at least six columns.")
  if(nrow(gtData) == 0) stop("gtData must have at least one row.")
  dimnames(gtData) <- NULL
  if(any(!is.finite(gtData))) stop("All values in gtData must be finite.")
  z.obs <- gtData[, 1]
  if(any(z.obs != 0 & z.obs != 1)){
    stop("Pooled test outcomes in column 1 of gtData must be either 0 or 1.")
  }
  psz <- gtData[, 2]
  if(any(psz < 1) || any((psz%%1) != 0)){
    stop("Pool sizes in column 2 of gtData must be positive integers.")
  }

  # This block tracks individuals assigned to each pool
  Memb <- gtData[ ,-(1:5), drop=FALSE]
  if(any(rowSums(Memb > 0) != psz)){
    stop("Pool sizes in column 2 must equal the number of positive individual IDs in each row of gtData.")
  }
  dup.fn <- function(x){
    x <- x[x > 0]
    anyDuplicated(x) > 0
  }
  dup.member <- apply(Memb, 1, dup.fn)
  if(any(dup.member)){
    stop("An individual ID cannot appear more than once in the same pool.")
  }
  all.ids <- Memb[Memb > 0]
  if(length(all.ids) == 0) stop("gtData must have at least one positive individual ID.")
  if(any((all.ids%%1)!= 0)) stop("Individual IDs in gtData must be integers.")
  ids <- sort(unique(all.ids))
  if(!all(ids == seq_len(max(ids)))) {
    stop("Individual IDs in gtData must be consecutive integers starting at 1.")
  }
  N <- max(ids)
  maxAssign <- max(as.numeric(table(all.ids)))
  ytm <- matrix(-9L,N,maxAssign)
  vec <- 1:nrow(gtData)
  for(d in 1:N){
    tid <- Memb==d
    store <- NULL
    for(i in 1:ncol(Memb)) {
      store <- c(store, vec[tid[ ,i]])
    }
    ytm[d,1:length(store)] <- sort(store)
  }

  ## Validating the EM arguments
  if(!is.numeric(nburn) || length(nburn) != 1L || !is.finite(nburn) || nburn < 0 || nburn %% 1 != 0){
    stop("nburn must be a non-negative integer.")
  }
  if(!is.numeric(ngit) || length(ngit) != 1L || !is.finite(ngit) || ngit <= 0 || ngit %% 1 != 0){
    stop("ngit must be a positive integer.")
  }
  if(!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) || maxit <= 0 || maxit %% 1 != 0){
    stop("maxit must be a positive integer.")
  }
  if(!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0){
    stop("tol must be a positive finite numeric value.")
  }
  if(!is.logical(tracing) || length(tracing) != 1L || is.na(tracing)){
    stop("tracing must be TRUE or FALSE.")
  }

  # Validating covariance & conf.level
  if(!is.logical(covariance) || length(covariance) != 1L || is.na(covariance)) stop("covariance must be TRUE or FALSE.")
  if(!is.numeric(conf.level) || length(conf.level) != 1L || !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1){
    stop("conf.level must be a number between 0 and 1.")
  }

  ## Generating individual true disease statuses
  ## at the initial parameter value p0
  Yt <- stats::rbinom(N,1,p0)
  Ytmat <- cbind(Yt,rowSums(ytm>0),ytm)

  ## Some global variables
  Ycol <- ncol(Ytmat)
  SeSp <- gtData[ ,3:4]
  if(any(SeSp < 0) || any(SeSp > 1)){
    stop("Sensitivity and specificity values in columns 3 and 4 of gtData must be between 0 and 1.")
  }
  Z <- gtData[ ,-(3:5)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  
  ## The total number of Gibbs iterates
  GI <- ngit + nburn
  
  ## Initial value of the parameter
  p1 <- p0
  p0 <- p0 + 2*tol
  s <- 0
  convergence <- 0
  
  ## EM algorithm starts here
  while(abs(p1-p0) > tol){ 
    s <- s + 1
    p0 <- p1
    U <- matrix(stats::runif(N*GI),nrow=N,ncol=GI)

    ## Gibbs sampling in Fortran to approximate the E-step
    res <- .Call("gbsonedhom_c",as.double(p0),as.integer(Ytmat),
                 as.integer(Z),as.integer(N),as.double(SeSp),as.integer(Ycol),
                 as.integer(Zrow),as.integer(Zcol),as.double(U),as.integer(GI),
                 as.integer(nburn), PACKAGE="groupTesting")
    temp <- sum( res )/ngit

    ## M-step: The parameter p is updated here
    p1 <- temp/N

    ## Trace the progress
    if(tracing){
      cat(s, p1, "\n")
    }
    ## Check convergence
    if(abs(p1-p0) <= tol){
      break
    }
    ## Terminate the EM algorithm
    ## if it exceeds max iteration
    if(s >= maxit){
      convergence <- 1
      break
    }
  }

  ## Covariance matrix in Fortran using Louis's (1982) method
  covr2 <- NULL
  if(covariance){
    U <- matrix(stats::runif(N*GI),nrow=N,ncol=GI)
    Info <- .Call("cvondknachom_c",as.double(p1),as.integer(Ytmat),
                  as.integer(Z),as.integer(N),as.double(SeSp),as.integer(Ycol),
                  as.integer(Zrow),as.integer(Zcol),as.double(U),as.integer(GI),
                  as.integer(nburn), PACKAGE="groupTesting")
    covr2 <- solve( Info )
  }

  ## Univariate Wald Inference
  pHat <- p1
  if(covariance){
    se <- sqrt(covr2)
    alternative <- "two.sided"
    ## Calculate the test statistic:
    z <- stats::qnorm(ifelse(alternative=="two.sided",
                     (1+conf.level)/2, conf.level))
    ## Find the confidence interval:
    CI <- c(pHat-z*se, pHat+z*se)
    ## Organizing and reporting the output:
    res <- data.frame(round(pHat, 3), round(se, 3), 
                      round(CI[1],3), round(CI[2],3) )
  }else{
    res <- data.frame(round(pHat, 3), NA, NA, NA) 
  }

  rownames(res) <- colnames(res) <- NULL
  rownames(res) <- "prop"
  colnames(res) <- c("Estimate", "Std.Err",                    
                     paste(conf.level*100, "%lower", sep=""),
                     paste(conf.level*100, "%upper", sep=""))

  # Output
  list("param"       = p1,
       "covariance"  = covr2,
	   "iterUsed"    = s,
       "convergence" = convergence,
	   "summary"     = res
	   )
}
