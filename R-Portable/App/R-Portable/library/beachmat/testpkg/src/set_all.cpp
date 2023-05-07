#include "set_all.h"
#include "get_all.h"

extern "C" {

// Set all rows.

SEXP set_row_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_row_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(optr.get(), order));
    END_RCPP
}

SEXP set_row_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);
    
    set_row_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_row_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(optr.get(), order));
    END_RCPP
}

SEXP set_row_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_row_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(optr.get(), order));
    END_RCPP
}

SEXP set_row_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::CharacterVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_row_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(optr.get(), order));
    END_RCPP
}

// Set all columns.

SEXP set_col_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    
    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_col_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(optr.get(), order));
    END_RCPP
}

SEXP set_col_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_col_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(optr.get(), order));
    END_RCPP
}

SEXP set_col_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_col_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(optr.get(), order));
    END_RCPP
}

SEXP set_col_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::CharacterVector>(ptr.get(), optr.get(), order);
    return Rcpp::List::create(optr->yield(), get_col_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(optr.get(), order));
    END_RCPP
}

// Set individual entries.

SEXP set_single_all_numeric (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);
    
    set_single_all<Rcpp::NumericVector>(ptr.get(), optr.get(), rorder, corder);
    return Rcpp::List::create(optr->yield(), get_single_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(optr.get(), rorder, corder));
    END_RCPP
}

SEXP set_single_all_integer (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);
    
    set_single_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), rorder, corder);
    return Rcpp::List::create(optr->yield(), get_single_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(optr.get(), rorder, corder));
    END_RCPP
}

SEXP set_single_all_logical (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_single_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), rorder, corder);
    return Rcpp::List::create(optr->yield(), get_single_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(optr.get(), rorder, corder));
    END_RCPP
}

SEXP set_single_all_character (SEXP in, SEXP rorder, SEXP corder) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_character_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_single_all<Rcpp::CharacterVector>(ptr.get(), optr.get(), rorder, corder);
    return Rcpp::List::create(optr->yield(), get_single_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(optr.get(), rorder, corder));
    END_RCPP
}

}
