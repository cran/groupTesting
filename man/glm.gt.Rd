% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regEM.R
\name{glm.gt}
\alias{glm.gt}
\title{EM Algorithm for Fitting Regression Models to Group Testing Data}
\usage{
glm.gt(
  beta0,
  gtData,
  X,
  g,
  dg = NULL,
  d2g = NULL,
  grdMethod = c("central", "forward", "backward"),
  covariance = FALSE,
  nburn = 2000,
  ngit = 5000,
  maxit = 200,
  tol = 0.001,
  tracing = TRUE,
  conf.level = 0.95,
  ...
)
}
\arguments{
\item{beta0}{An initial value for the regression coefficients.}

\item{gtData}{A matrix or data.frame consisting of the pooled test outcomes and other information from a group testing application. Needs to be specified as shown in the example below.}

\item{X}{The design matrix.}

\item{g}{An inverse link function in the GLM family.}

\item{dg}{The first derivate of \code{g}. When NULL, a finite-difference approximation will be used.}

\item{d2g}{The second derivate of \code{g}. When NULL, a finite-difference approximation will be used.}

\item{grdMethod}{The finite-difference approximation method to be used for \code{dg} and \code{d2g}. See 'Details'.}

\item{covariance}{When TRUE, the covariance matrix is calculated at the MLE.}

\item{nburn}{The number of initial Gibbs iterates to be discarded.}

\item{ngit}{The number of Gibbs iterates to be used in the E-step after discarding \code{nburn} iterates as a burn-in period.}

\item{maxit}{The maximum number of EM steps (iterations) allowed in the EM algorithm.}

\item{tol}{Convergence tolerance used in the EM algorithm.}

\item{tracing}{When TRUE, progress in the EM algorithm is displayed.}

\item{conf.level}{Confidence level to be used for the Wald confidence interval.}

\item{...}{Further arguments to be passed to \code{\link{optim}}. See 'Details'.}
}
\value{
A list with components:
\item{param}{The MLE of the regression coefficients.}
\item{covariance}{Estimated covariance matrix for the regression coefficients.}
\item{iterUsed}{The number of EM iterations needed for convergence.}
\item{convergence}{0 if the EM algorithm converges successfully and 1 if the iteration limit \code{maxit} has been reached.}
\item{summary}{Estimation summary with Wald confidence interval.}
}
\description{
This function implements an expectation-maximization (EM) algorithm to fit regression models to group testing data, where pooled responses are related to individual covariates through a link function in the generalized linear model (GLM) family. The EM algorithm, which is outlined in Warasi (2023), finds the maximum likelihood estimate (MLE) for the vector of regression coefficients, \strong{beta}. The EM algorithm can model pooling data observed from \strong{any} group testing protocol used in practice, including hierarchical and array testing (Kim et al., 2007).
}
\details{
\code{gtData} must be specified as follows. Columns 1-5 consist of the pooled test outcomes (0 for negative and 1 for positive), pool sizes, pool-specific sensitivities, pool-specific specificities, and assay identification (ID) numbers, respectively. From column 6 onward, the pool member ID numbers need to be specified. Note that the ID numbers must start with 1 and increase consecutively up to \code{N}, the total number of individuals tested. \strong{For smaller pools, incomplete ID numbers must be filled out by -9 or any non-positive numbers} as shown in the example below. The design matrix \code{X} consists of invidual covariate information, such as age, sex, and symptoms, of the pool members located in column 6 onward.\tabular{cccccccccccc}{
    \tab Z \tab psz \tab Se \tab Sp \tab Assay \tab Mem1 \tab Mem2 \tab Mem3 \tab Mem4 \tab Mem5 \tab Mem6 \cr
   Pool:1 \tab 1 \tab 6 \tab 0.90 \tab 0.92 \tab 1 \tab 1 \tab 2 \tab 3 \tab 4 \tab 5 \tab 6 \cr
   Pool:2 \tab 0 \tab 6 \tab 0.90 \tab 0.92 \tab 1 \tab 7 \tab 8 \tab 9 \tab 10 \tab 11 \tab 12 \cr
   Pool:3 \tab 1 \tab 2 \tab 0.95 \tab 0.96 \tab 2 \tab 1 \tab 2 \tab -9 \tab -9 \tab -9 \tab -9 \cr
   Pool:4 \tab 0 \tab 2 \tab 0.95 \tab 0.96 \tab 2 \tab 3 \tab 4 \tab -9 \tab -9 \tab -9 \tab -9 \cr
   Pool:5 \tab 1 \tab 2 \tab 0.95 \tab 0.96 \tab 2 \tab 5 \tab 6 \tab -9 \tab -9 \tab -9 \tab -9 \cr
   Pool:6 \tab 0 \tab 1 \tab 0.92 \tab 0.90 \tab 3 \tab 1 \tab -9 \tab -9 \tab -9 \tab -9 \tab -9 \cr
   Pool:7 \tab 1 \tab 1 \tab 0.92 \tab 0.90 \tab 3 \tab 2 \tab -9 \tab -9 \tab -9 \tab -9 \tab -9 \cr
   Pool:8 \tab 0 \tab 1 \tab 0.92 \tab 0.90 \tab 3 \tab 5 \tab -9 \tab -9 \tab -9 \tab -9 \tab -9 \cr
   Pool:9 \tab 0 \tab 1 \tab 0.92 \tab 0.90 \tab 3 \tab 6 \tab -9 \tab -9 \tab -9 \tab -9 \tab -9 \cr
}


This is an example of \code{gtData}, where 12 individuals are assigned to 2 non-overlapping initial pools and then tested based on the 3-stage hierarchical protocol. The test outcomes, \code{Z}, from 9 pools are in column 1. In three stages, different pool sizes (6, 2, and 1), sensitivities, specificities, and assays are used. The ID numbers of the pool members are shown in columns 6-11. The row names and column names are not required. Note that the EM algorithm can accommodate any group testing data including those described in Kim et al. (2007). For individual testing data, the pool size in column 2 is 1 for all pools.

\code{X} is an \eqn{N}x\eqn{k} design matrix, where each column represents a vector of individual covariate values. For an intercept model, the first column values must be 1. The column (covariate) names of X, such as 'age' and 'sex', will be displayed in the estimation summary. When column names are missing (NULL), the names that will be displayed by default are 'Intercept', 'x1', 'x2', and so on.

The EM algorithm implements a Gibbs sampler to approximate the expectation in the E-step. Under each EM iteration, \code{ngit} Gibbs samples are retained for these purposes after discarding the initial \code{nburn} samples.

\code{g} relates the pooled responses Z (column 1 in \code{gtData}) to \code{X}. \code{dg} and \code{d2g} can be specified analogously. These characteristics can be obtained from \code{\link{glmLink}} for the common links: logit, probit, and complementary log-log.

\code{grdMethod} is used only when dg and d2g are NULL, where a finite-difference approximation is implemented by the function \code{fderiv} from the package 'pracma'.

The optimization routine \code{\link{optim}} is used to complete the M-step with the default method 'Nelder-Mead'. The argument ... allows the user to change the default method as well as other arguments in \code{\link{optim}}.

The covariance matrix is calculated by an appeal to the missing data principle and the method outlined in Louis (1982).
}
\examples{

library(groupTesting)

## To illustrate 'glm.gt', we use data simulated  
## by the functions 'hier.gt.simulation' and 'array.gt.simulation'.

## Note: The simulated data-structures are consistent  
## with the data-structure required for 'gtData'.

## Example 1: MLE from 3-stage hierarchical group testing data.
## The data used is simulated by 'hier.gt.simulation'. 

N <- 200              # Sample size
S <- 3                # 3-stage hierarchical testing
psz <- c(6,2,1)       # Pool sizes used in stages 1-3
Se <- c(.95,.95,.98)  # Sensitivities in stages 1-3
Sp <- c(.95,.98,.96)  # Specificities in stages 1-3
assayID <- c(1,2,3)   # Assays used in stages 1-3
param.t <- c(-3,2,1)  # The TRUE parameter to be estimated

# Simulating covariates:
set.seed(123)
x1 <- rnorm(N, mean=0, sd=0.75)
x2 <- rbinom(N, size=1, prob=0.5)
X <- cbind(1, x1, x2)
colnames( X ) <- c("Intercept", "Predictor 1", "Predictor 2")
# Note: Because the 1st column of X is 1, intercept model will be fit.

# Specifying logit inverse link:
g <- function(t){exp(t)/(1+exp(t))}  
pReg <- g(X\%*\%param.t)

# Simulating test responses:
gtOut <- hier.gt.simulation(N,pReg,S,psz,Se,Sp,assayID)$gtData

# Fitting the model (with intercept):
param0 <- param.t + 0.2   # Initial value
res <- glm.gt(beta0=param0,gtData=gtOut,X=X,
              g=g,dg=NULL,d2g=NULL,
              grdMethod="central",covariance=TRUE,
              nburn=2000,ngit=5000,maxit=200,
              tol=1e-03,tracing=TRUE,conf.level=0.95)

# Note: Because dg and d2g are NULL (i.e., the exact derivatives
#       are not given), numerical derivatives are used.

# Estimation results:
# > res

# $param
# [1] -2.840802  1.992916  0.677176

# $covariance
#            [,1]        [,2]        [,3]
# [1,]  0.2134439 -0.10147555 -0.16693776
# [2,] -0.1014756  0.16855122  0.02997113
# [3,] -0.1669378  0.02997113  0.26324589

# $iterUsed
# [1] 10

# $convergence
# [1] 0

# $summary
#             Estimate Std.Err 95\%lower 95\%upper
# Intercept     -2.841   0.462   -3.746   -1.935
# Predictor 1    1.993   0.411    1.188    2.798
# Predictor 2    0.677   0.513   -0.328    1.683

## Example 2: MLE from two-dimensional array testing data.
## The data used is simulated by 'array.gt.simulation'. 

N <- 200            # Sample size
protocol <- "A2"     # 2-stage array without testing initial master pool
n <- 5               # Row/column size
Se <- c(0.95, 0.95)  # Sensitivities
Sp <- c(0.98, 0.98)  # Specificities
assayID <- c(1, 1)   # The same assay in both stages
param <- c(-4,1,1)   # The TRUE parameter to be estimated

# Simulating data:
set.seed(123)
x1 <- runif(N)
x2 <- rnorm(N, mean=0, sd=0.5)
x3 <- rbinom(N, size=1, prob=0.5)
X <- cbind(x1, x2, x3)
# Note: Because the 1st column of X is not 1, 
#       the model without intercept will be fit.

# Finding g, dg, and d2g from the function 'glmLink':  
res0 <- glmLink(fn.name="logit")  
g <- res0$g            # Logit inverse link g()
dg <- res0$dg          # The exact first derivate of g
d2g <- res0$d2g        # The exact second derivate of g
pReg <- g(X\%*\%param)   # Individual probabilities
gtOut <- array.gt.simulation(N,pReg,protocol,n,Se,Sp,assayID)$gtData

# Fitting the model (without intercept):
param0 <- param + 0.2 
res <- glm.gt(beta0=param0,gtData=gtOut,X=X,g=g,                  
              dg=dg,d2g=d2g,covariance=TRUE,
              nburn=2000,ngit=5000,maxit=200,
              tol=1e-03,tracing=TRUE,conf.level=0.95)
print(res)

\donttest{
## Example 3: MLE from non-overlapping initial pooled responses.
## The data used is simulated by 'hier.gt.simulation'.

## Note: With initial pooled responses, our MLE is equivalent  
## to the MLE in Vansteelandt et al. (2000).

N <- 1000             # Sample size
psz <- 5              # Pool size
S <- 1                # 1-stage testing
Se <- 0.95            # Sensitivity
Sp <- 0.99            # Specificity
assayID <- 1          # Assay used for all pools
param <- c(-3,2,1)    # The TRUE parameter to be estimated

# Simulating data:
set.seed(123)
x1 <- rnorm(N, mean=0, sd=0.75)
x2 <- rbinom(N, size=1, prob=0.5)
X <- cbind(1, x1, x2)

# Finding g, dg, and d2g by the function 'glmLink':  
res0 <- glmLink(fn.name="probit")  # Probit link
g <- res0$g
dg <- res0$dg
d2g <- res0$d2g
pReg <- g(X\%*\%param)
gtOut <- hier.gt.simulation(N,pReg,S,psz,Se,Sp,assayID)$gtData

# Fitting the model:
param0 <- param + 0.2
res <- glm.gt(beta0=param0,gtData=gtOut,X=X,g=g,                 
                 dg=dg,d2g=d2g,covariance=TRUE,
                 nburn=2000,ngit=5000,maxit=200,
                 tol=1e-03,tracing=TRUE,conf.level=0.95)
print(res)

## Example 4: MLE from individual (one-by-one) testing data.
## The data used is simulated by 'hier.gt.simulation'.

N <- 1000             # Sample size
psz <- 1              # Pool size 1 (i.e., individual testing)
S <- 1                # 1-stage testing
Se <- 0.95            # Sensitivity
Sp <- 0.99            # Specificity
assayID <- 1          # Assay used for all pools
param <- c(-3,2,1)    # The TRUE parameter to be estimated

# Simulating data:
set.seed(123)
x1 <- rnorm(N, mean=0, sd=0.75)
x2 <- rbinom(N, size=1, prob=0.5)
X <- cbind(1, x1, x2)
g <- function(t){exp(t)/(1+exp(t))}  # Inverse logit 
pReg <- g(X\%*\%param)
gtOut <- hier.gt.simulation(N,pReg,S,psz,Se,Sp,assayID)$gtData

# Fitting the model:
param0 <- param + 0.2 
res <- glm.gt(beta0=param0,gtData=gtOut,
              X=X,g=g,dg=NULL,d2g=NULL,
              grdMethod="central",covariance=TRUE,
              nburn=2000,ngit=5000,maxit=200,
              tol=1e-03,tracing=TRUE,conf.level=0.95)
print(res)

## Example 5: Using pooled testing data.

# Pooled test outcomes:
Z <- c(1, 0, 1, 0, 1, 0, 1, 0, 0)     

# Design matrix, X:
x1 <- c(0.8,1.2,0.4,1.5,1.8,1.8,0.1,1.6,0.2,0.2,1.8,0.2)
x2 <- c(31,56,45,64,26,47,22,60,35,41,32,41)
X <- cbind(x1, x2)

# Pool sizes used:
psz <- c(6, 6, 2, 2, 2, 1, 1, 1, 1)

# Pool-specific Se & Sp:
Se <- c(.90, .90, .95, .95, .95, .92, .92, .92, .92)
Sp <- c(.92, .92, .96, .96, .96, .90, .90, .90, .90)

# Assays used:
Assay <- c(1, 1, 2, 2, 2, 3, 3, 3, 3)

# Pool members:
Memb <- rbind( 
   c(1, 2,  3,  4,  5,  6),
   c(7, 8,  9, 10, 11, 12),
   c(1, 2, -9, -9, -9, -9),
   c(3, 4, -9, -9, -9, -9),
   c(5, 6, -9, -9, -9, -9),
   c(1,-9, -9, -9, -9, -9),
   c(2,-9, -9, -9, -9, -9),
   c(5,-9, -9, -9, -9, -9),
   c(6,-9, -9, -9, -9, -9)
)
# The data-structure suited for 'gtData':
gtOut <- cbind(Z, psz, Se, Sp, Assay, Memb)

# Fitting the model with logit link:
g <- function(t){exp(t)/(1+exp(t))} 
param0 <- c(0, 0)      
res <- glm.gt(beta0=param0,gtData=gtOut,X=X,                  
              g=g,dg=NULL,d2g=NULL,
              grdMethod="central",covariance=TRUE,
              nburn=2000,ngit=5000,maxit=200,
              tol=1e-03,tracing=TRUE,conf.level=0.95)
print(res)
}

}
\references{
Kim HY, Hudgens M, Dreyfuss J, Westreich D, and Pilcher C. (2007). Comparison of Group Testing Algorithms for Case Identification in the Presence of Testing Error. \emph{Biometrics}, 63:1152-1163.

Louis T. (1982). Finding the Observed Information Matrix when Using the EM algorithm. \emph{Journal of the Royal Statistical Society: Series B}, 44:226-233.

Vansteelandt S, Goetghebeur E, and Verstraeten T. (2000). Regression Models for Disease Prevalence with Diagnostic Tests on Pools of Serum Samples. \emph{Biometrics}, 56:1126-1133.

Warasi M. (2023). groupTesting: An R Package for Group Testing Estimation. \emph{Communications in Statistics-Simulation and Computation}, 52:6210-6224.
}
\seealso{
\code{\link{hier.gt.simulation}} and \code{\link{array.gt.simulation}} for group testing data simulation, and \code{\link{prop.gt}} for estimation of a disease prevalence from group testing data.
}
