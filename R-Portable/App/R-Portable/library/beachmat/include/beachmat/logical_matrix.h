#ifndef BEACHMAT_LOGICAL_MATRIX_H
#define BEACHMAT_LOGICAL_MATRIX_H

#include "input/LIN_matrix.h"
#include "output/LIN_output.h"
#include "utils/dispatch.h"

#include <memory>
#include <string>

namespace beachmat {

/*********
 * INPUT *
 *********/

/* Virtual base class for logical matrices. */

typedef lin_matrix<int, Rcpp::LogicalVector> logical_matrix;

std::unique_ptr<logical_matrix> create_logical_matrix_internal(const Rcpp::RObject&, bool); 

/* Simple logical matrix */

typedef simple_lin_matrix<int, Rcpp::LogicalVector> simple_logical_matrix;

/* lgeMatrix */

template<>
inline std::string dense_reader<int, Rcpp::LogicalVector>::get_class() { return "lgeMatrix"; }

typedef dense_lin_matrix<int, Rcpp::LogicalVector> dense_logical_matrix;

/* lgCMatrix */

template<>
inline int Csparse_reader<int, Rcpp::LogicalVector>::get_empty() { return 0; }

template<>
inline std::string Csparse_reader<int, Rcpp::LogicalVector>::get_class() { return "lgCMatrix"; }

typedef Csparse_lin_matrix<int, Rcpp::LogicalVector> Csparse_logical_matrix;

/* DelayedMatrix */

template<>
inline std::unique_ptr<logical_matrix> delayed_lin_reader<int, Rcpp::LogicalVector>::generate_seed(Rcpp::RObject incoming) {
    return create_logical_matrix_internal(incoming, false);
}

typedef delayed_lin_matrix<int, Rcpp::LogicalVector> delayed_logical_matrix;

/* Unknown matrix */

typedef unknown_lin_matrix<int, Rcpp::LogicalVector> unknown_logical_matrix;

/* External matrix */

template<>
inline std::string external_reader_base<int, Rcpp::LogicalVector>::get_type() { return "logical"; } 

typedef external_lin_matrix<int, Rcpp::LogicalVector> external_logical_matrix;

/* Dispatcher */

inline std::unique_ptr<logical_matrix> create_logical_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) {
        std::string ctype=get_class_name(incoming);
        if (ctype=="lgeMatrix") { 
            return std::unique_ptr<logical_matrix>(new dense_logical_matrix(incoming));
        } else if (ctype=="lgCMatrix") { 
            return std::unique_ptr<logical_matrix>(new Csparse_logical_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<logical_matrix>(new delayed_logical_matrix(incoming));
        } else if (has_external_support("logical", incoming)) {
            return std::unique_ptr<logical_matrix>(new external_logical_matrix(incoming));
        }
        return std::unique_ptr<logical_matrix>(new unknown_logical_matrix(incoming));
    } 
    quit_on_df(incoming);
    return std::unique_ptr<logical_matrix>(new simple_logical_matrix(incoming));
}

inline std::unique_ptr<logical_matrix> create_logical_matrix(const Rcpp::RObject& incoming) {
    return create_logical_matrix_internal(incoming, true);
}

template<>
inline std::unique_ptr<logical_matrix> create_matrix<logical_matrix>(const Rcpp::RObject& incoming) {
    return create_logical_matrix(incoming);
}

/**********
 * OUTPUT *
 **********/

/* Virtual base class for output logical matrices. */

typedef lin_output<int, Rcpp::LogicalVector> logical_output;

/* Simple output logical matrix */

typedef simple_lin_output<int, Rcpp::LogicalVector> simple_logical_output;

/* Sparse output logical matrix */

template<>
inline int Csparse_writer<int, Rcpp::LogicalVector>::get_empty() { return 0; }

template<>
inline std::string Csparse_writer<int, Rcpp::LogicalVector>::get_class() { return "lgCMatrix"; }

typedef sparse_lin_output<int, Rcpp::LogicalVector> sparse_logical_output;

/* Simple output logical matrix */

template<>
inline std::string external_writer_base<int, Rcpp::LogicalVector>::get_type() { return "logical"; } 

typedef external_lin_output<int, Rcpp::LogicalVector> external_logical_output;

/* Output dispatchers */

inline std::unique_ptr<logical_output> create_logical_output(int nrow, int ncol, const output_param& param) {
    auto pkg=param.get_package();

    if (pkg=="Matrix") {
        if (param.get_class()=="lgCMatrix") {
            return std::unique_ptr<logical_output>(new sparse_logical_output(nrow, ncol));
        }
    } else if (param.is_external_available("logical")) {
        return std::unique_ptr<logical_output>(new external_logical_output(nrow, ncol, pkg, param.get_class()));
    }

    return std::unique_ptr<logical_output>(new simple_logical_output(nrow, ncol));
}

template<>
inline std::unique_ptr<logical_output> create_output<logical_output>(int nrow, int ncol, const output_param& param) {
    return create_logical_output(nrow, ncol, param);
}

}

#endif
