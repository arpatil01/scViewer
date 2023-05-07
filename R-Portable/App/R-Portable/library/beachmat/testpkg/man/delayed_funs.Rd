\name{Delayed functions}
\alias{delayed_funs}

\title{Create DelayedMatrices}
\description{Create instances of \linkS4class{DelayedMatrix} objects with different sets of delayed operations.}

\usage{
delayed_funs(basefun, DELAYED_FUN)
}

\arguments{
\item{basefun}{Function that creates a non-DelayedMatrix matrix representation.}
\item{DELAYED_FUN}{Function that alters the values of the individual elements in a matrix representation.}
}

\value{
A list of functions where each function generates a DelayedMatrix instance with a different type of delayed operation.
}

\details{
DelayedMatrix objects where the only delayed operations are subsetting or transposition will be handled specially by \pkg{beachmat}.
This necessitates checking of different types of DelayedMatrix subsetting inputs to ensure that this handling is correct.
Any other delayed operations will cause the object to be handled by block processing.
}

\author{
Aaron Lun
}

\examples{
library(Matrix)
FUNS <- delayed_funs(rsparsematrix, DELAYED_FUN=function(x) x * 2)
FUNS[[1]](10, 20, 0.2)
}
