#include "set_slice.h"
#include "get_slice.h"

extern "C" {

// Set all rows.

SEXP set_row_slice_numeric (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_slice<Rcpp::NumericVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_row_slice_integer (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_slice<Rcpp::IntegerVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_row_slice_logical (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_slice<Rcpp::LogicalVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_row_slice_character (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_slice<Rcpp::CharacterVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(optr.get(), order, subset));
    END_RCPP
}

// Set all columns.

SEXP set_col_slice_numeric (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_slice<Rcpp::NumericVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_col_slice_integer (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_slice<Rcpp::IntegerVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_col_slice_logical (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_slice<Rcpp::LogicalVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_col_slice_character (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_slice<Rcpp::CharacterVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(optr.get(), order, subset));
    END_RCPP
}

}
