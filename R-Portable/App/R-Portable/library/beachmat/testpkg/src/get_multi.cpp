#include "get_multi.h"

extern "C" {

// Get all rows.

SEXP get_multirow_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_multirow_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_multirow_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_multirow_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_multirow_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_multirow_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_multirow_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_multirow_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order);
    END_RCPP
}

// Get all columns.

SEXP get_multicol_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_multicol_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_multicol_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_multicol_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_multicol_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_multicol_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_multicol_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_multicol_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order);
    END_RCPP
}

// Get row slice.

SEXP get_multirow_slice_numeric (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_multirow_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order, bound);
    END_RCPP
}

SEXP get_multirow_slice_integer (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_multirow_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order, bound);
    END_RCPP
}

SEXP get_multirow_slice_logical (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_multirow_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order, bound);
    END_RCPP
}

SEXP get_multirow_slice_character (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_multirow_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order, bound);
    END_RCPP
}

// Get column slice.

SEXP get_multicol_slice_numeric (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_multicol_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order, bound);
    END_RCPP
}

SEXP get_multicol_slice_integer (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_multicol_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order, bound);
    END_RCPP
}

SEXP get_multicol_slice_logical (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_multicol_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order, bound);
    END_RCPP
}

SEXP get_multicol_slice_character (SEXP in, SEXP order, SEXP bound) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_multicol_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order, bound);
    END_RCPP
}

}
