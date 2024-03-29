% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtcode.R
\name{array.gt.simulation}
\alias{array.gt.simulation}
\title{Simulating Array-Based Group Testing Data}
\usage{
array.gt.simulation(
  N,
  p = 0.1,
  protocol = c("A2", "A2M"),
  n,
  Se,
  Sp,
  assayID,
  Yt = NULL
)
}
\arguments{
\item{N}{The number of individuals to be tested.}

\item{p}{A vector of length N consisting of individual disease probabilities.}

\item{protocol}{Either "A2" or "A2M", where "A2" ("A2M") refers to the two-dimensional array without (with) testing the members of an array as a single pooled sample.}

\item{n}{The row (or column) size of the arrays.}

\item{Se}{A vector of assay sensitivities.}

\item{Sp}{A vector of assay specificities.}

\item{assayID}{A vector of assay identification numbers.}

\item{Yt}{A vector of individual true disease statuses.}
}
\value{
A list with components:
\item{gtData}{The simulated group testing data.}
\item{testsExp}{The number of tests expended in the simulation.}
}
\description{
This function simulates two-dimensional array-based group testing data.
}
\details{
We consider the array testing protocol outlined in Kim et al. (2007). Under this protocol, \eqn{N} individuals are assigned to \eqn{m} non-overlapping \eqn{n}-by-\eqn{n} matrices such that \eqn{N=mn^2}. From each matrix, \eqn{n} pools are formed using the row specimens and another \eqn{n} pools are formed using the column specimens. In stage 1, the \eqn{2n} pools are tested. In stage 2, individual testing is used for case identification according to the strategy described in Kim et al. (2007). This is a 2-stage protocol called \emph{Square Array without Master Pool Testing} and denoted by \eqn{A2(n:1)} in Kim et al. (2007). A variant (3-stage protocol) is also presented in Kim et al. (2007) which employs testing the \eqn{n^2} array members together as an initial pooled unit before implementing the 2-stage array. If the initial pooled test is negative, the procedure stops (i.e., the 2-stage array is not needed). However, if the pooled test is positive, the 2-stage protocol is used as before. This 3-stage approach is called \emph{Square Array with Master Pool Testing} and is denoted by \eqn{A2(n^2:n:1)}. See Kim et al. (2007) for more details.

\code{N} should be divisible by the array size \eqn{n^2}. When not divisible, the remainder individuals are tested one by one (i.e., individual testing).

\code{p} is a vector of individual disease probabilities. When all individuals have the same probability of disease, say, 0.10, \code{p} can be specified as rep(0.10, N) or p=0.10.

For "A2" and "A2M", the pool sizes used are \code{c(n, 1)} and \code{c(n^2, n, 1)}, respectively.

For "A2", \code{Se} is \code{c(Se1, Se2)}, where \code{Se1} is the sensitivity of the assay used for both row and column pools, and \code{Se2} is the sensitivity of the assay used for individual testing. For "A2M", \code{Se} is \code{c(Se1, Se2, Se3)}, where \code{Se1} is for the initial array pool, \code{Se2} is for the row and column pools, and \code{Se3} is for individual testing. \code{Sp} is specified in the same manner.

For "A2", \code{assayID} is \code{c(1, 1)} when the same assay is used for row/column pool testing as well as for individual testing, and assayID is \code{c(1, 2)} when assay 1 is used for row/column pool testing and assay 2 is used for individual testing. In the same manner, \code{assayID} is specified for "A2M" as \code{c(1, 1, 1)}, \code{c(1, 2, 3)}, and in many other ways.

When available, the individual true disease statuses (1 for positive and 0 for negative) can be used in simulating the group testing data through argument \code{Yt}. When an input is entered for \code{Yt}, argument \code{p} will be ignored.
}
\examples{

library(groupTesting)

## Example 1: Square Array without Master Pool Testing (i.e., 2-Stage Array)
N <- 48              # Sample size
protocol <- "A2"     # 2-stage array
n <- 4               # Row/column size
Se <- c(0.95, 0.95)  # Sensitivities in stages 1-2
Sp <- c(0.98, 0.98)  # Specificities in stages 1-2
assayID <- c(1, 1)   # The same assay in both stages

# (a) Homogeneous population
pHom <- 0.10         # Overall prevalence
array.gt.simulation(N=N,p=pHom,protocol=protocol,n=n,Se=Se,Sp=Sp,assayID=assayID)

# Alternatively, the individual true statuses can be used as: 
yt <- rbinom( N, size=1, prob=0.1 )
array.gt.simulation(N=N,protocol=protocol,n=n,Se=Se,Sp=Sp,assayID=assayID,Yt=yt)

# (b) Heterogeneous population (regression)
param <- c(-3,2,1)
x1 <- rnorm(N, mean=0, sd=.75)
x2 <- rbinom(N, size=1, prob=0.5)
X <- cbind(1, x1, x2)
pReg <- exp(X\%*\%param)/(1+exp(X\%*\%param)) # Logit
array.gt.simulation(N=N,p=pReg,protocol=protocol,n=n,Se=Se,Sp=Sp,assayID=assayID)

# The above examples with different assays
Se <- c(0.95, 0.98)
Sp <- c(0.97, 0.99)
assayID <- c(1, 2)
array.gt.simulation(N,pHom,protocol,n,Se,Sp,assayID)
array.gt.simulation(N,pReg,protocol,n,Se,Sp,assayID)

## Example 2: Square Array with Master Pool Testing (i.e., 3-Stage Array)
N <- 48
protocol <- "A2M"
n <- 4
Se <- c(0.95, 0.95, 0.95)
Sp <- c(0.98, 0.98, 0.98)
assayID <- c(1, 1, 1)    # The same assay in 3 stages

# (a) Homogeneous population
pHom <- 0.10
array.gt.simulation(N,pHom,protocol,n,Se,Sp,assayID)

# (b) Heterogeneous population (regression)
param <- c(-3,2,1)
x1 <- rnorm(N, mean=0, sd=.75)
x2 <- rbinom(N, size=1, prob=0.5)
X <- cbind(1, x1, x2)
pReg <- exp(X\%*\%param)/(1+exp(X\%*\%param)) # Logit
array.gt.simulation(N,pReg,protocol,n,Se,Sp,assayID)

# The above examples with different assays:
Se <- c(0.95, 0.98, 0.98)
Sp <- c(0.97, 0.98, 0.92)
assayID <- 1:3
array.gt.simulation(N,pHom,protocol,n,Se,Sp,assayID)
array.gt.simulation(N,pReg,protocol,n,Se,Sp,assayID)

}
\references{
Kim HY, Hudgens M, Dreyfuss J, Westreich D, and Pilcher C (2007). Comparison of Group Testing Algorithms for Case Identification in the Presence of Testing Error. \emph{Biometrics}, 63(4), 1152–1163.
}
\seealso{
\code{\link{hier.gt.simulation}} for simulation of the hierarchical group testing data.
}
