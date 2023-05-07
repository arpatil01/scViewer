#include "get_all.h"

extern "C" {

// Get all rows.

SEXP get_row_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_row_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_row_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_row_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_row_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_row_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_row_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_row_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order);
    END_RCPP
}

// Get all columns.

SEXP get_col_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_col_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_col_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_col_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_col_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_col_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_col_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_col_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order);
    END_RCPP
}

// Get individual entries.

SEXP get_single_all_numeric (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_single_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), rorder, corder);
    END_RCPP
}

SEXP get_single_all_integer (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_single_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), rorder, corder);
    END_RCPP
}

SEXP get_single_all_logical (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_single_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), rorder, corder);
    END_RCPP
}

SEXP get_single_all_character (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_single_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), rorder, corder);
    END_RCPP
}

}
