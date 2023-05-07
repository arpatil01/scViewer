#ifndef RATICATE_CSPARSEMATRIX_HPP
#define RATICATE_CSPARSEMATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include "tatami/ext/ArrayView.hpp"

namespace raticate { 

template<typename Data = double, typename Index = int, class V>
Parsed<Data, Index> parse_CSparseMatrix(Rcpp::RObject seed, const V& val) {
    auto dims = parse_dims(seed.slot("Dim"));

    Rcpp::IntegerVector p(seed.slot("p"));
    tatami::ArrayView<int> pview(static_cast<const int*>(p.begin()), p.size());
    Rcpp::IntegerVector i(seed.slot("i"));
    tatami::ArrayView<int> iview(static_cast<const int*>(i.begin()), i.size());

    typedef typename std::remove_const<typename std::remove_reference<decltype(val[0])>::type>::type Value;
    tatami::ArrayView<Value> vview(static_cast<const Value*>(val.begin()), val.size());

    Parsed<Data, Index> output;
    output.matrix.reset(
        new tatami::CompressedSparseMatrix<false, Data, Index, decltype(vview), decltype(iview), decltype(pview)>(
            dims.first, 
            dims.second, 
            std::move(vview), 
            std::move(iview), 
            std::move(pview)
        )
    );
    output.contents = Rcpp::List::create(i, val, p); // protect views from GC of the underlying memory.

    return output;
}

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_dgCMatrix(Rcpp::RObject seed) {
    Rcpp::NumericVector y(seed.slot("x"));
    return parse_CSparseMatrix(seed, y);
}

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_lgCMatrix(Rcpp::RObject seed) {
    Rcpp::LogicalVector y(seed.slot("x"));
    return parse_CSparseMatrix(seed, y);
}

}

#endif
