#ifndef RATICATE_SIMPLEMATRIX_HPP
#define RATICATE_SIMPLEMATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include "tatami/ext/ArrayView.hpp"

namespace raticate { 

template<typename Data = double, typename Index = int, class V>
Parsed<Data, Index> parse_simple_matrix_internal(const V& y) {
    Parsed<Data, Index> output;

    typedef typename std::remove_const<typename std::remove_reference<decltype(y[0])>::type>::type Value;
    tatami::ArrayView view(static_cast<const Value*>(y.begin()), y.size());
    output.matrix.reset(new tatami::DenseColumnMatrix<double, int, decltype(view)>(y.rows(), y.cols(), std::move(view)));

    output.contents = Rcpp::List::create(y);
    return output;
}

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_simple_matrix(const Rcpp::RObject& seed) {
    Parsed<Data, Index> output;

    if (seed.sexp_type() == REALSXP) {
        Rcpp::NumericMatrix y(seed);
        output = parse_simple_matrix_internal(y);
    } else if (seed.sexp_type() == INTSXP) {
        Rcpp::IntegerMatrix y(seed);
        output = parse_simple_matrix_internal(y);
    } else if (seed.sexp_type() == LGLSXP) {
        Rcpp::LogicalMatrix y(seed);
        output = parse_simple_matrix_internal(y);
    }

    return output;
}

}

#endif
