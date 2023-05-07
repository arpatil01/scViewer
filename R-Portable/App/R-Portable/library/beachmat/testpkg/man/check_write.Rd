\name{Check writes}
\alias{check_write_all}
\alias{check_write_slice}
\alias{check_write_varslice}
\alias{check_write_type}
\alias{check_write_class}
\alias{check_write_indexed}
\alias{check_write_errors}

\title{Check writing with \pkg{beachmat}}
\description{Check that \pkg{beachmat} can successfully write to a matrix representation.}

\usage{
check_write_all(FUN, ..., mode, out.class=NULL)

check_write_slice(FUN, ..., mode, out.class=NULL)

check_write_varslice(FUN, ..., mode, out.class=NULL)

check_write_type(FUN, ..., mode, out.class=NULL)

check_write_class(test.mat, expected)

check_write_indexed(FUN, ..., mode, out.class=NULL)

check_write_errors(FUN, ..., mode, out.class=NULL)
}

\arguments{
    \item{FUN}{A function that generates an instance of a matrix representation.}
    \item{...}{Further arguments to pass to \code{FUN}.}
    \item{mode}{String specifying whether the matrix is \code{"integer"}, \code{"logical"}, \code{"numeric"} or \code{"character"}.}
    \item{out.class}{String specifying the expected class of the output matrix.
        If \code{NULL}, it defaults to the class of the input matrix representation.}
    \item{test.mat}{An instance of a matrix representation.}
    \item{expected}{String specifying the class of \code{test.mat}, as detected by \pkg{beachmat}.}
}

\value{
If writing is correct, \code{NULL} is invisibly returned.
Otherwise, any errors in writeing are handled via the \pkg{testthat} framework. 
}

\details{
These functions will check whether the matrix representation can be correctly written by \pkg{beachmat} in C++ code.
This is done by making a copy of an input matrix and checking that the output instance has the same representation and values.
The rows and columns are shuffled during access to avoid trivial equality by simply returning the input.

Assuming that we have a matrix \code{x} containing \dQuote{default} values for each data type (e.g., \code{0} or \code{""}):
\itemize{
\item \code{check_write_all} will check writing to each row (\code{x[i,] <-} for scalar \code{i}), column (\code{x[,j] <-} for scalar \code{j}) and element \code{x[i,j] <-}).
\item \code{check_write_slice} will check writing to a subset of each row (\code{x[i,a:b] <-}) or column (\code{x[a:b,j] <-}).
\item \code{check_write_varslice} is the same as \code{check_write_slice} except that a different subset is used for each row/column.
\item \code{check_write_type} will check the integer/double conversion upon single row/column writing in integer, logical or numeric matrices.
Nothing is done if \code{x} is a character matrix.
\item \code{check_write_class} will check that the class of \code{x} is correctly identified.
\item \code{check_write_indexed} will check for indexed writing to each column (\code{x[i,v] <-} for vector \code{v}) or row (\code{x[v,j] <-}).
\item \code{check_write_errors} will check that error conditions are correctly raised for invalid writing requests to \code{x}.
}

These functions will also simultaneously check matrix access from the same C++ classes used to write to the output matrix.
}

\author{
Aaron Lun
}

\examples{
FUN <- function(nr=10, nc=20) matrix(rnorm(nr*nc), nr, nc)
check_write_all(FUN=FUN, mode="numeric")
}
