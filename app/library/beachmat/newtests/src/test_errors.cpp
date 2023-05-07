#include "beachmat3/beachmat.h"
#include <map>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject test_sparse_writer1(int type) {
    std::map<std::pair<int, int>, double> thingy;

    if (type==1) {
        thingy[std::make_pair(-1, 0)] = 1;
    } else if (type==2) {
        thingy[std::make_pair(0, -1)] = 1;
    } else if (type==3) {
        thingy[std::make_pair(0, 100)] = 1;
    } else if (type==4) {
        thingy[std::make_pair(100, 0)] = 1;
    }

    return beachmat::as_gCMatrix<Rcpp::NumericVector>(1, 1, thingy);
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject test_sparse_writer2(Rcpp::RObject mat, Rcpp::NumericVector replacement) {
    return beachmat::as_gCMatrix<Rcpp::NumericVector>(mat, replacement);
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject test_sparse_writer3() {
    std::deque<std::pair<std::pair<int, int>, double> > thingy;
    thingy.push_back(std::make_pair(std::make_pair(1, 1), 1));
    thingy.push_back(std::make_pair(std::make_pair(0, 0), 1));
    return beachmat::as_gCMatrix<Rcpp::NumericVector>(2, 2, thingy);
}
