#include "get_const.h"

extern "C" {

// Get all const columns.

SEXP get_const_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::create_numeric_matrix(in);
    return get_const_all<Rcpp::NumericMatrix>(mat.get(), order);
    END_RCPP
}

SEXP get_const_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::create_integer_matrix(in);
    return get_const_all<Rcpp::IntegerMatrix>(mat.get(), order);
    END_RCPP
}

SEXP get_const_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::create_logical_matrix(in);
    return get_const_all<Rcpp::LogicalMatrix>(mat.get(), order);
    END_RCPP
}

SEXP get_const_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::create_character_matrix(in);
    return get_const_all<Rcpp::CharacterMatrix>(mat.get(), order);
    END_RCPP
}

// Get const column slices.

SEXP get_const_slice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_numeric_matrix(in);
    return get_const_slice<Rcpp::NumericMatrix>(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_const_slice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_integer_matrix(in);
    return get_const_slice<Rcpp::IntegerMatrix>(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_const_slice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_logical_matrix(in);
    return get_const_slice<Rcpp::LogicalMatrix>(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_const_slice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_character_matrix(in);
    return get_const_slice<Rcpp::CharacterMatrix>(mat.get(), order, bounds);
    END_RCPP
}

// Get variable const column slices.

SEXP get_const_varslice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_numeric_matrix(in);
    return get_const_varslice(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_const_varslice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_integer_matrix(in);
    return get_const_varslice(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_const_varslice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_logical_matrix(in);
    return get_const_varslice(mat.get(), order, bounds);
    END_RCPP
}

SEXP get_const_varslice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::create_character_matrix(in);
    return get_const_varslice(mat.get(), order, bounds);
    END_RCPP
}

}
