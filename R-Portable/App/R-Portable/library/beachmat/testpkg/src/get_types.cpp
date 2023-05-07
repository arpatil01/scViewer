#include "get_all.h"

extern "C" {

// Type check.

SEXP get_type(SEXP in) {
    BEGIN_RCPP
    std::string out=beachmat::translate_type(beachmat::find_sexp_type(in));
    return Rf_mkString(out.c_str());
    END_RCPP
}

// Row conversion.

SEXP get_row_numeric_to_logical (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_row_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_row_numeric_to_integer (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_row_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_row_logical_to_numeric (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_row_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_row_logical_to_integer (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_row_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_row_integer_to_numeric (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_row_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_row_integer_to_logical (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_row_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), ordering);
    END_RCPP
}

// Column conversion.

SEXP get_col_numeric_to_logical (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_col_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_col_numeric_to_integer (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_col_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_col_logical_to_numeric (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_col_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_col_logical_to_integer (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_col_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_col_integer_to_numeric (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_col_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), ordering);
    END_RCPP
}

SEXP get_col_integer_to_logical (SEXP in, SEXP ordering) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_col_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), ordering);
    END_RCPP
}

}
