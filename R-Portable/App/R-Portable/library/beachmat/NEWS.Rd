\name{beachmatnews}
\title{beachmat News}
\encoding{UTF-8}

\section{Version 2.14.0}{\itemize{
\item Added the \pkg{tatami} C++ library for LinkingTo from other packages.
}}

\section{Version 2.8.0}{\itemize{
\item Improve the efficiency of sparse row subsetting in non-DelayedArray \code{rowBlockApply()}.

\item Avoid overhead of \pkg{DelayedArray} block processing when DelayedArray is pristine and the type is supported.

\item Migrated \code{whichNonZero()} from \pkg{scuttle}.

\item Added \code{toCsparse()} to make it easier to convert SparseArraySeeds to CsparseMatrixes.

\item Added \code{realizeFileBackedMatrix()} to, well, realize a DelayedMatrix with file-backed components.
}}

\section{Version 2.6.0}{\itemize{
\item Added version 3 of the C++ API.
This relies on \pkg{DelayedArray} block processing to handle the various abstractions,
and provides convenient C++ abstractions for the different input types.
More methods for non-zero extraction are provided for sparse algorithms,
and the class hierarchy handles flexibility of the input types rather than requiring user-side templating.

\item Added the \code{colBlockApply()} and \code{rowBlockApply()} functions for more convenient and efficient column/row-block processing.

\item Added doxygen documentation for the C++ API.
}}

\section{Version 2.0.0}{\itemize{
\item Converted beachmat into a header-only library to simplify linking, avoid ABI incompatibilities.

\item Minor API change for external native support to avoid using C++ classes with C-style linkage.
}}

\section{Version 1.4.0}{\itemize{
\item Removed native support for RleMatrix and packed symmetric matrices.

\item Added multi-row/column getters.

\item Added mechanism for native support of arbitrary developer-defined matrices.

\item Switched to row/colGrid() for defining chunks in unsupported matrices.
}}

\section{Version 1.2.0}{\itemize{
\item Changed environment variable to BEACHMAT_RPATH for consistency with other packages.

\item Added native support for transposition and subsetting in DelayedMatrix objects.

\item Added support for chunk-by-chunk realization of otherwise unsupported matrices, including DelayedMatrix objects with other delayed operations.

\item Added the get_const_col_indexed() method for input matrices, especially fast for sparse representations. 

\item Added the set_col_indexed() and set_row_indexed() methods for output matrices.

\item Updated vignettes and expanded the test suite.
}}

\section{Version 1.0.0}{\itemize{
\item
New package beachmat, which provides a C++ API for handling various R matrix types.
}}

