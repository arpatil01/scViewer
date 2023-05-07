#include "get_slice.h"

extern "C" {

// Get row slice.

SEXP get_row_slice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_row_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_row_slice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_row_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_row_slice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_row_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_row_slice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_row_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

// Get column slice.

SEXP get_col_slice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_col_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_col_slice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_col_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_col_slice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_col_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_col_slice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_col_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

}
