#include "beachtest.h"

extern "C" {

SEXP get_class_integer(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(incoming);
    return Rcpp::StringVector::create(ptr->get_class());	
    END_RCPP
}

SEXP get_class_numeric(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(incoming);
    return Rcpp::StringVector::create(ptr->get_class());	
    END_RCPP
}

SEXP get_class_logical(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(incoming);
    return Rcpp::StringVector::create(ptr->get_class());	
    END_RCPP
}

SEXP get_class_character(SEXP incoming) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(incoming);
    return Rcpp::StringVector::create(ptr->get_class());	
    END_RCPP
}

}

