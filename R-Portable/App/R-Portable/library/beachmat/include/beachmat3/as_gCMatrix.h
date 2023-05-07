#ifndef BEACHMAT_AS_GCMATRIX_H
#define BEACHMAT_AS_GCMATRIX_H

/**
 * @file as_gCMatrix.h
 *
 * Functions to create `*gCMatrix` instances.
 */

#include "Rcpp.h"
#include <map>

namespace beachmat {

/** 
 * Generate a `*gCMatrix` of the appropriate class for the specified `Rcpp::Vector` class in `V`.
 *
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @tparam V An `Rcpp::Vector` class.
 * Currently only `Rcpp::NumericVector` and `Rcpp::LogicalVector` are supported.
 *
 * @return An `Rcpp::S4` object of the appropriate `*gCMatrix` class for `V`,
 * either `dgCMatrix` or `lgCMatrix` for double-precision and logical data, respectively.
 * All other choices of `V` will trigger a compile-time error.
 */
template <class V>
inline Rcpp::S4 generate_gCMatrix () { // could also use = delete here.
    static_assert(sizeof(V)==0, "unsupported specialization of generate_gCMatrix");
    return Rcpp::S4("lgCMatrix");
}

template <>
inline Rcpp::S4 generate_gCMatrix<Rcpp::LogicalVector> () {
    return Rcpp::S4("lgCMatrix");
}

template <>
inline Rcpp::S4 generate_gCMatrix<Rcpp::NumericVector> () {
    return Rcpp::S4("dgCMatrix");
}

/**
 * Create a `*gCMatrix` from a triplet-formatted store of non-zero entries.
 * Best used when the number of non-zero entries is not known in advance.
 *
 * @tparam V An `Rcpp::Vector` class, to be used as the `x` slot in the output `*gCMatrix`.
 * Only `Rcpp::NumericVector` and `Rcpp::LogicalVector` are supported.
 * @tparam S An appropriate triplet store, see comments below.
 *
 * @param nr Number of rows.
 * @param nc Number of columns.
 * @param store A sorted triplet store.
 *
 * @details
 * For an element `x` in the triplet store, we should obtain the zero-based _column_ index from `x.first.first`;
 * the zero-based _row_ index from `x.first.second`; and the value from `x.second`.
 * The store itself can be any container that supports the usual STL `.begin()`, `.end()` and `.size()` operations.
 * Elements in the store should be sorted by the column index and row index.
 * (If multiple elements have the same row/column indices, the last element in the store is used.)
 * One can naturally obtain an appropriate store with a `std::map<std::pair<int, int>, double>`, which sorts for us as well.
 * Alternatively, we can create a store manually with `std::deque<std::pair<std::pair<int, int>, double> >`.
 *
 * @return A `dgCMatrix` or `lgCMatrix` instance (depending on `V`) containing all entries in `store`.
 */
template <class V, class S>
inline Rcpp::RObject as_gCMatrix (int nr, int nc, const S& store) {
    auto mat = generate_gCMatrix<V>();
    mat.slot("Dim") = Rcpp::IntegerVector::create(nr, nc);

    const size_t nnzero = store.size();
    Rcpp::IntegerVector i(nnzero);
    V x(nnzero);
    Rcpp::IntegerVector p(nc + 1, 0);

    auto xIt=x.begin();
    auto iIt=i.begin();
    auto sIt = store.begin();

    int counter = 0;
    int lastcol = 0, lastrow = 0;

    for (int c = 1; c <= nc; ++c) {
        while (sIt != store.end() && (sIt->first).first < c) {
            const auto& curcol = (sIt->first).first;
            const auto& currow = (sIt->first).second;

            if (currow >= nr || currow < 0) {
                throw std::runtime_error("entries in 'store' refer to out-of-range rows");
            }
            if (curcol < 0) {
                throw std::runtime_error("entries in 'store' refer to out-of-range columns");
            }
            if (lastcol > curcol || (lastcol == curcol && lastrow > currow)) {
                throw std::runtime_error("entries in 'store' are not sorted");
            } 

            (*xIt) = (sIt->second);
            (*iIt) = currow;
            lastcol = curcol;
            lastrow = currow;

            ++xIt;
            ++iIt;
            ++sIt;
            ++counter;
        }
        p[c] = counter;
    }
    
    if (static_cast<size_t>(counter) != store.size()) {
        throw std::runtime_error("entries in 'store' refer to out-of-range columns");
    }

    mat.slot("p")=p;
    mat.slot("i")=i;
    mat.slot("x")=x;

    return SEXP(mat);
}

/**
 * Create a `*gCMatrix` by modifying the non-zero values (but not their order or position) in an existing `*gCMatrix`.
 *
 * @tparam V An `Rcpp::Vector` class, to be used as the `x` slot in the output `*gCMatrix`.
 * Only `Rcpp::NumericVector` and `Rcpp::LogicalVector` are supported.
 * @tparam T Type of data in the triplet store.
 * It should be possible to cast this to the storage type of `V`.
 *
 * @param old A `*gCMatrix` object.
 * This does not need to be of a type corresponding to `V`.
 * @param x A vector of non-zero values to use to replace the existing `x` slot in `old`.
 * The row index and column identity oeach value is assumed to be the same as the non-zero values already in `old`.
 *
 * @return A `dgCMatrix` or `lgCMatrix` instance (depending on `V`) 
 * containing indexing information from `old` but updated with the non-zero values in `x`.
 */
template <class V>
inline Rcpp::RObject as_gCMatrix (Rcpp::RObject old, V x) {
    // Needs to be regenerated from scratch in case of a type conversion.
    auto mat = generate_gCMatrix<V>();

    mat.slot("Dim") = old.slot("Dim");
    mat.slot("p") = old.slot("p");
    Rcpp::IntegerVector oldi = old.slot("i");
    mat.slot("i") = oldi;

    if (x.size() != oldi.size()) {
        throw std::runtime_error("inconsistent number of non-zero entries");
    }
    mat.slot("x") = x;


    return SEXP(mat);
}

}

#endif
