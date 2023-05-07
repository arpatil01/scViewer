#include "get_errors.h"

extern "C" { 

SEXP get_errors_integer (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    get_errors<Rcpp::IntegerVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_errors_logical (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    get_errors<Rcpp::LogicalVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_errors_numeric (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    get_errors<Rcpp::NumericVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_errors_character (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    get_errors<Rcpp::StringVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_multi_errors_integer (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    get_multi_errors<Rcpp::IntegerVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_multi_errors_logical (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    get_multi_errors<Rcpp::LogicalVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_multi_errors_numeric (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    get_multi_errors<Rcpp::NumericVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}

SEXP get_multi_errors_character (SEXP in, SEXP mode) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    get_multi_errors<Rcpp::StringVector>(ptr.get(), mode);
    return Rf_ScalarLogical(1);
    END_RCPP
}
}
