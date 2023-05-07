\name{BiocNeighborsnews}
\title{BiocNeighbors News}
\encoding{UTF-8}

\section{Version 1.10.0}{\itemize{
\item Migrated \code{findMutualNN()} from \pkg{batchelor}.

\item Vendored the \pkg{RcppAnnoy} headers for greater reproducibility.

\item Added a \code{distance="Cosine"} option for all algorithms.
}}

\section{Version 1.4.0}{\itemize{
\item Allow memory-efficient retrieval of the distance to the furthest neighbors.

\item Added a warn.ties= argument to turn off tie-related warnings in the KMKNN and VP tree algorithms.

\item Return neighbor counts in rangeFind*() and rangeQuery*() functions when get.index=FALSE and get.distance=FALSE.
}}

\section{Version 1.2.0}{\itemize{
\item findNeighbors() and queryNeighbors() now accept a vector of point-specific thresholds.

\item Added an VP tree implemention in findVptree(), queryVptree(), buildVptree().
Supported dispatch to these methods in the generics.

\item Added a HNSW implementation in findHnsw(), queryHnsw(), buildHnsw().
Supported dispatch to these methods in the generics.

\item Renamed buildNNIndex() to buildIndex().

\item Converted findNeighbors() and queryNeighbors() into S4 methods.
Created specific rangeFind*() and rangeQuery*() functions for KMKNN and VP tree algorithms.

\item Modified AnnoyIndex class to hold the original data matrix.
Created bnorder(), bndata() generics to obtain (possibly reordered) input matrices for all indexes.

\item Supported Manhattan distance searching in all algorithms.
}}

\section{Version 1.0.0}{\itemize{
\item
New package BiocNeighbors, for k-means-based k-nearest neighbor detection.
}}

