#ifndef BEACHMAT_INTEGER_MATRIX_H
#define BEACHMAT_INTEGER_MATRIX_H

#include "input/LIN_matrix.h"
#include "output/LIN_output.h"
#include "utils/dispatch.h"

#include <memory>
#include <string>

namespace beachmat {

/*************
 *** INPUT ***
 *************/

/* Virtual base class for integer matrices. */

typedef lin_matrix<int, Rcpp::IntegerVector> integer_matrix;

std::unique_ptr<integer_matrix> create_integer_matrix_internal(const Rcpp::RObject&, bool);

/* Simple integer matrix */

typedef simple_lin_matrix<int, Rcpp::IntegerVector> simple_integer_matrix;

/* DelayedMatrix */

template<>
inline std::unique_ptr<integer_matrix> delayed_lin_reader<int, Rcpp::IntegerVector>::generate_seed(Rcpp::RObject incoming) {
    return create_integer_matrix_internal(incoming, false);
}

typedef delayed_lin_matrix<int, Rcpp::IntegerVector> delayed_integer_matrix;

/* Unknown matrix */

typedef unknown_lin_matrix<int, Rcpp::IntegerVector> unknown_integer_matrix;

/* External matrix */

template<>
inline std::string external_reader_base<int, Rcpp::IntegerVector>::get_type() { return "integer"; }

typedef external_lin_matrix<int, Rcpp::IntegerVector> external_integer_matrix;

/* Dispatcher */

inline std::unique_ptr<integer_matrix> create_integer_matrix_internal(const Rcpp::RObject& incoming, bool delayed) {
    if (incoming.isS4()) { 
        std::string ctype=get_class_name(incoming);
        if (delayed && ctype=="DelayedMatrix") {
            return std::unique_ptr<integer_matrix>(new delayed_integer_matrix(incoming));
        } else if (has_external_support("integer", incoming)) {
            return std::unique_ptr<integer_matrix>(new external_integer_matrix(incoming));
        }
        return std::unique_ptr<integer_matrix>(new unknown_integer_matrix(incoming));
    }
    quit_on_df(incoming);
    return std::unique_ptr<integer_matrix>(new simple_integer_matrix(incoming));
}

inline std::unique_ptr<integer_matrix> create_integer_matrix(const Rcpp::RObject& incoming) {
    return create_integer_matrix_internal(incoming, true);
}

template<>
inline std::unique_ptr<integer_matrix> create_matrix<integer_matrix>(const Rcpp::RObject& incoming) {
    return create_integer_matrix(incoming);
}

/**************
 *** OUTPUT ***
 **************/

/* Virtual base class for output integer matrices. */

typedef lin_output<int, Rcpp::IntegerVector> integer_output;

/* Simple output integer matrix */

typedef simple_lin_output<int, Rcpp::IntegerVector> simple_integer_output;

/* External output integer matrix */

template<>
inline std::string external_writer_base<int, Rcpp::IntegerVector>::get_type() { return "integer"; }

typedef external_lin_output<int, Rcpp::IntegerVector> external_integer_output;

/* Output dispatchers */

inline std::unique_ptr<integer_output> create_integer_output(int nrow, int ncol, const output_param& param) {
    if (param.is_external_available("integer")) { 
        return std::unique_ptr<integer_output>(new external_integer_output(nrow, ncol, param.get_package(), param.get_class()));
    }
     
    return std::unique_ptr<integer_output>(new simple_integer_output(nrow, ncol));
}

template<>
inline std::unique_ptr<integer_output> create_output<integer_output>(int nrow, int ncol, const output_param& param) {
    return create_integer_output(nrow, ncol, param);
}

}

#endif
