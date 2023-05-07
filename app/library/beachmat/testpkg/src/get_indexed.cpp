#include "get_indexed.h"

extern "C" {

// Get all const columns.

SEXP get_indexed_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::create_numeric_matrix(in);
    return get_indexed_all<Rcpp::NumericMatrix>(mat.get(), order);
    END_RCPP
}

SEXP get_indexed_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::create_logical_matrix(in);
    return get_indexed_all<Rcpp::LogicalMatrix>(mat.get(), order);
    END_RCPP
}

// Get const column slices.

SEXP get_indexed_slice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_numeric_matrix(in);
    return get_indexed_slice<Rcpp::NumericMatrix>(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_indexed_slice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_logical_matrix(in);
    return get_indexed_slice<Rcpp::LogicalMatrix>(mat.get(), order, bounds);
    END_RCPP
}

// Get variable const column slices.

SEXP get_indexed_varslice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_numeric_matrix(in);
    return get_indexed_varslice(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_indexed_varslice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_logical_matrix(in);
    return get_indexed_varslice(mat.get(), order, bounds);
    END_RCPP
}

}
