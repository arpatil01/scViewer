\name{Check reads}
\alias{check_read_all}
\alias{check_read_slice}
\alias{check_read_varslice}
\alias{check_read_multi}
\alias{check_read_type}
\alias{check_read_class}
\alias{check_read_const}
\alias{check_read_indexed}
\alias{check_read_errors}

\title{Check reading with \pkg{beachmat}}
\description{Check that \pkg{beachmat} can successfully read a matrix representation.}

\usage{
check_read_all(FUN, ..., mode)

check_read_slice(FUN, ..., mode)

check_read_varslice(FUN, ..., mode)

check_read_multi(FUN, ..., mode)

check_read_type(FUN, ..., mode)

check_read_class(test.mat, mode, expected)

check_read_const(FUN, ..., mode)

check_read_indexed(FUN, ..., mode)

check_read_errors(FUN, ..., mode)
}

\arguments{
    \item{FUN}{A function that generates an instance of a matrix representation.}
    \item{...}{Further arguments to pass to \code{FUN}.}
    \item{mode}{String specifying whether the matrix is \code{"integer"}, \code{"logical"}, \code{"numeric"} or \code{"character"}.}
    \item{test.mat}{An instance of a matrix representation.}
    \item{expected}{String specifying the class of \code{test.mat}, as detected by \pkg{beachmat}.}
}

\value{
If reading is correct, \code{NULL} is invisibly returned.
Otherwise, any errors in reading are handled via the \pkg{testthat} framework. 
}

\details{
These functions will check whether the matrix representation can be correctly accessed by \pkg{beachmat} in C++ code.
For a given matrix \code{x}:
\itemize{
\item \code{check_read_all} will check access to each row (\code{x[i,]} for scalar \code{i}), column (\code{x[,j]} for scalar \code{j}) and element \code{x[i,j]}).
\item \code{check_read_slice} will check access to a subset of each row (\code{x[i,a:b]}) or column (\code{x[a:b,j]}).
\item \code{check_read_varslice} is the same as \code{check_read_slice} except that a different subset is used for each row/column.
\item \code{check_read_multi} will check access to multiple rows (\code{x[i1:i2,]}) or columns (\code{x[,j1:j2]}).
It will also check access to subsets of multiple rows (e.g., \code{x[i1:i2, a:b]}) and columns.
\item \code{check_read_type} will check the integer/double conversion upon single row/column access in integer, logical or numeric matrices.
Nothing is done if \code{x} is a character matrix.
\item \code{check_read_class} will check that the class of \code{x} is correctly identified.
\item \code{check_read_const} will check for constant copy-free access to columns of an ordinary matrix.
For all other matrices, it checks for with-copy column access.
\item \code{check_read_indexed} will check for constant copy-free access to row indices and values for columns of a \linkS4class{dgCMatrix}.
For all other matrices, it checks for with-copy column access.
\item \code{check_read_errors} will check that error conditions are correctly raised for invalid access requests to \code{x}.
}
}

\author{
Aaron Lun
}

\examples{
FUN <- function(nr=10, nc=20) matrix(rnorm(nr*nc), nr, nc)
check_read_all(FUN=FUN, mode="numeric")
}
