#include "beachtest.h"

extern "C" {

SEXP set_class_by_sexp(SEXP incoming) {
    BEGIN_RCPP
    beachmat::output_param op(incoming);
    return Rcpp::StringVector::create(op.get_class());
    END_RCPP
}

}

