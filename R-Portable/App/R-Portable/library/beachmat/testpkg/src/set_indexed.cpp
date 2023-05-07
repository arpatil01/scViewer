#include "set_indexed.h"

extern "C" {

// Set row indices.

SEXP set_row_indexed_integer(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    return set_row_indexed<Rcpp::IntegerVector>(optr.get(), index1, index2);
    END_RCPP
}

SEXP set_row_indexed_logical(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    return set_row_indexed<Rcpp::LogicalVector>(optr.get(), index1, index2);
    END_RCPP
}

SEXP set_row_indexed_numeric(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    
    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    return set_row_indexed<Rcpp::NumericVector>(optr.get(), index1, index2);
    END_RCPP
}

SEXP set_row_indexed_character(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    return set_row_indexed<Rcpp::StringVector>(optr.get(), index1, index2);
    END_RCPP
}

// Set col indices.

SEXP set_col_indexed_integer(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);
    return set_col_indexed<Rcpp::IntegerVector>(optr.get(), index1, index2);
    END_RCPP
}

SEXP set_col_indexed_logical(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);
    return set_col_indexed<Rcpp::LogicalVector>(optr.get(), index1, index2);
    END_RCPP
}

SEXP set_col_indexed_numeric(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);
    return set_col_indexed<Rcpp::NumericVector>(optr.get(), index1, index2);
    END_RCPP
}

SEXP set_col_indexed_character(SEXP in, SEXP index1, SEXP index2) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    return set_col_indexed<Rcpp::StringVector>(optr.get(), index1, index2);
    END_RCPP
}

}
