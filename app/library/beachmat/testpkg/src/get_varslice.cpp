#include "get_varslice.h"

extern "C" {

// Get row slice.

SEXP get_row_varslice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_row_varslice<Rcpp::NumericVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_row_varslice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_row_varslice<Rcpp::IntegerVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_row_varslice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_row_varslice<Rcpp::LogicalVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_row_varslice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_row_varslice<Rcpp::CharacterVector>(ptr.get(), order, bounds);
    END_RCPP
}

// Get column slice.

SEXP get_col_varslice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_col_varslice<Rcpp::NumericVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_col_varslice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_col_varslice<Rcpp::IntegerVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_col_varslice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_col_varslice<Rcpp::LogicalVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_col_varslice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_col_varslice<Rcpp::CharacterVector>(ptr.get(), order, bounds);
    END_RCPP
}

}
