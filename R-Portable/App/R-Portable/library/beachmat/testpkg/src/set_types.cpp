#include "get_all.h"
#include "set_all.h"

extern "C" {

// Row conversion.

SEXP set_row_numeric_via_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_numeric_via_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_integer_via_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_integer_via_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_logical_via_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_row_logical_via_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_row_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

// Column conversion.

SEXP set_col_numeric_via_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_numeric_via_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_numeric_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_integer_via_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::LogicalVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_integer_via_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_integer_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_logical_via_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::IntegerVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

SEXP set_col_logical_via_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);

    beachmat::output_param op(ptr->get_class(), ptr->get_package());
    auto optr=beachmat::create_logical_output(ptr->get_nrow(), ptr->get_ncol(), op);

    set_col_all<Rcpp::NumericVector>(ptr.get(), optr.get(), order);
    return optr->yield();
    END_RCPP
}

}
