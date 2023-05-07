#ifndef BEACHMAT_NUMERIC_MATRIX_H
#define BEACHMAT_NUMERIC_MATRIX_H

#include "input/LIN_matrix.h"
#include "output/LIN_output.h"
#include "utils/dispatch.h"

#include <memory>
#include <string>

namespace beachmat { 

/*********
 * INPUT *
 *********/ 

/* Virtual base class for numeric matrices. */

typedef lin_matrix<double, Rcpp::NumericVector> numeric_matrix;

std::unique_ptr<numeric_matrix> create_numeric_matrix_internal(const Rcpp::RObject&, bool); 

/* Simple numeric matrix */

typedef simple_lin_matrix<double, Rcpp::NumericVector> simple_numeric_matrix;

/* dgeMatrix */

template<>
inline std::string dense_reader<double, Rcpp::NumericVector>::get_class() { return "dgeMatrix"; }

typedef dense_lin_matrix<double, Rcpp::NumericVector> dense_numeric_matrix;

/* dgCMatrix */

template<>
inline double Csparse_reader<double, Rcpp::NumericVector>::get_empty() { return 0; }

template<>
inline std::string Csparse_reader<double, Rcpp::NumericVector>::get_class() { return "dgCMatrix"; }

typedef Csparse_lin_matrix<double, Rcpp::NumericVector> Csparse_numeric_matrix;

/* DelayedMatrix */

template<>
inline std::unique_ptr<numeric_matrix> delayed_lin_reader<double, Rcpp::NumericVector>::generate_seed(Rcpp::RObject incoming) {
    return create_numeric_matrix_internal(incoming, false);
}

typedef delayed_lin_matrix<double, Rcpp::NumericVector> delayed_numeric_matrix;

/* Unknown matrix */

typedef unknown_lin_matrix<double, Rcpp::NumericVector> unknown_numeric_matrix;

/* External matrix */

template<>
inline std::string external_reader_base<double, Rcpp::NumericVector>::get_type() { return "numeric"; }

typedef external_lin_matrix<double, Rcpp::NumericVector> external_numeric_matrix;

/* Dispatcher */

inline std::unique_ptr<numeric_matrix> create_numeric_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) {
        std::string ctype=get_class_name(incoming);
        if (ctype=="dgeMatrix") { 
            return std::unique_ptr<numeric_matrix>(new dense_numeric_matrix(incoming));
        } else if (ctype=="dgCMatrix") { 
            return std::unique_ptr<numeric_matrix>(new Csparse_numeric_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<numeric_matrix>(new delayed_numeric_matrix(incoming));
        } else if (has_external_support("numeric", incoming)) {
            return std::unique_ptr<numeric_matrix>(new external_numeric_matrix(incoming));
        }
        return std::unique_ptr<numeric_matrix>(new unknown_numeric_matrix(incoming));
    } 
    quit_on_df(incoming);
    return std::unique_ptr<numeric_matrix>(new simple_numeric_matrix(incoming));
}

inline std::unique_ptr<numeric_matrix> create_numeric_matrix(const Rcpp::RObject& incoming) { 
    return create_numeric_matrix_internal(incoming, true);
}

template<>
inline std::unique_ptr<numeric_matrix> create_matrix<numeric_matrix>(const Rcpp::RObject& incoming) {
    return create_numeric_matrix(incoming);
}

/**********
 * OUTPUT *
 **********/ 

/* Virtual base class for output numeric matrices. */

typedef lin_output<double, Rcpp::NumericVector> numeric_output;

/* Simple output numeric matrix */

typedef simple_lin_output<double, Rcpp::NumericVector> simple_numeric_output;

/* Sparse output numeric matrix */

template<>
inline double Csparse_writer<double, Rcpp::NumericVector>::get_empty() { return 0; }

template<>
inline std::string Csparse_writer<double, Rcpp::NumericVector>::get_class() { return "dgCMatrix"; }

typedef sparse_lin_output<double, Rcpp::NumericVector> sparse_numeric_output;

/* External output numeric matrix */

template<>
inline std::string external_writer_base<double, Rcpp::NumericVector>::get_type() { return "numeric"; }

typedef external_lin_output<double, Rcpp::NumericVector> external_numeric_output;

/* Output dispatchers */

inline std::unique_ptr<numeric_output> create_numeric_output(int nrow, int ncol, const output_param& param) {
    auto pkg=param.get_package();

    if (pkg=="Matrix") {
        if (param.get_class()=="dgCMatrix") {
            return std::unique_ptr<numeric_output>(new sparse_numeric_output(nrow, ncol));
        }
    } else if (param.is_external_available("numeric")) {
        return std::unique_ptr<numeric_output>(new external_numeric_output(nrow, ncol, pkg, param.get_class()));
    }

    return std::unique_ptr<numeric_output>(new simple_numeric_output(nrow, ncol));
}

template<>
inline std::unique_ptr<numeric_output> create_output<numeric_output>(int nrow, int ncol, const output_param& param) {
    return create_numeric_output(nrow, ncol, param);
}

}

#endif
