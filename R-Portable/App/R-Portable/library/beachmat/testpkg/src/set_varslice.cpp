#include "set_varslice.h"
#include "get_varslice.h"

extern "C" {

// Set all rows.

SEXP set_row_varslice_numeric (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    
    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_varslice<Rcpp::NumericVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_varslice<Rcpp::NumericVector>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_row_varslice_integer (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_varslice<Rcpp::IntegerVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_varslice<Rcpp::IntegerVector>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_row_varslice_logical (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_varslice<Rcpp::LogicalVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_varslice<Rcpp::LogicalVector>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_row_varslice_character (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_varslice<Rcpp::CharacterVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_row_varslice<Rcpp::CharacterVector>(optr.get(), order, subset));
    END_RCPP
}

// Set all columns.

SEXP set_col_varslice_numeric (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_varslice<Rcpp::NumericVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_varslice<Rcpp::NumericVector>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_col_varslice_integer (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_varslice<Rcpp::IntegerVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_varslice<Rcpp::IntegerVector>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_col_varslice_logical (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_varslice<Rcpp::LogicalVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_varslice<Rcpp::LogicalVector>(optr.get(), order, subset));
    END_RCPP
}

SEXP set_col_varslice_character (SEXP in, SEXP order, SEXP subset) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_varslice<Rcpp::CharacterVector>(ptr.get(), optr.get(), order, subset);
    return Rcpp::List::create(optr->yield(), get_col_varslice<Rcpp::CharacterVector>(optr.get(), order, subset));
    END_RCPP
}

}
