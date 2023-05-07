#include "aaron_matrix.h"
#include "exports.h"

// Demonstrating with integer matrices.

typedef AaronMatrix<int, Rcpp::IntegerVector, Rcpp::IntegerMatrix> AaronIntMat;

// Constructor, destructors and clones.

void * AaronMatrix_integer_input_create (SEXP incoming) {
    return static_cast<void*>(new AaronIntMat(incoming));
}

void AaronMatrix_integer_input_destroy (void * ptr) {
    delete static_cast<AaronIntMat*>(ptr);
    return;
}

void * AaronMatrix_integer_input_clone (void * ptr) {
    AaronIntMat* old=static_cast<AaronIntMat*>(ptr);
    return static_cast<void*>(new AaronIntMat(*old));
}

// Basic getters

void AaronMatrix_integer_input_dim(void* ptr, size_t* nr, size_t* nc){ 
    AaronIntMat* thing=static_cast<AaronIntMat*>(ptr);
    *nr=thing->get_nrow();
    *nc=thing->get_ncol();
    return;
}

void AaronMatrix_integer_input_get(void * ptr, size_t r, size_t c, int* val) {
    *val=static_cast<AaronIntMat*>(ptr)->get(r, c);	
    return;
}

void AaronMatrix_integer_input_getRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_row(r, *out, first, last);
    return;
}

void AaronMatrix_integer_input_getCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_col(c, *out, first, last);
    return;
}

void AaronMatrix_integer_input_getRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_row(r, *out, first, last);
    return;
}

void AaronMatrix_integer_input_getCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Multi getters

void AaronMatrix_integer_input_getRows_integer(void * ptr, Rcpp::IntegerVector::iterator* r, size_t n, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_rows(*r, n, *out, first, last);
    return;
}

void AaronMatrix_integer_input_getCols_integer(void * ptr, Rcpp::IntegerVector::iterator* c, size_t n, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_cols(*c, n, *out, first, last);
    return;
}

void AaronMatrix_integer_input_getRows_numeric(void * ptr, Rcpp::IntegerVector::iterator* r, size_t n, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_rows(*r, n, *out, first, last);
    return;
}

void AaronMatrix_integer_input_getCols_numeric(void * ptr, Rcpp::IntegerVector::iterator* c, size_t n, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntMat*>(ptr)->get_cols(*c, n, *out, first, last);
    return;
}
