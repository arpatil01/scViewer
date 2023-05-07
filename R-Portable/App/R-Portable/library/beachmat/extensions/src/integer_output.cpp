#include "aaron_output.h"
#include "exports.h"

// Demonstrating with integer matrices.

typedef AaronOutput<int, Rcpp::IntegerVector, Rcpp::IntegerMatrix> AaronIntOut;

// Constructor, destructors and clones.

void * AaronMatrix_integer_output_create (size_t nr, size_t nc) {
    return static_cast<void*>(new AaronIntOut(nr, nc));
}

void AaronMatrix_integer_output_destroy (void * ptr) {
    delete static_cast<AaronIntOut*>(ptr);
    return;
}

void * AaronMatrix_integer_output_clone (void * ptr) {
    AaronIntOut* old=static_cast<AaronIntOut*>(ptr);
    return static_cast<void*>(new AaronIntOut(*old));
}

SEXP AaronMatrix_integer_output_yield (void * ptr) {
    return static_cast<AaronIntOut*>(ptr)->yield();
}

// Basic getters

void AaronMatrix_integer_output_get(void * ptr, size_t r, size_t c, int* val) {
    *val=static_cast<AaronIntOut*>(ptr)->get(r, c);	
    return;
}

void AaronMatrix_integer_output_getRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void AaronMatrix_integer_output_getCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

void AaronMatrix_integer_output_getRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void AaronMatrix_integer_output_getCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Basic setters

void AaronMatrix_integer_output_set(void * ptr, size_t r, size_t c, int* val) {
    static_cast<AaronIntOut*>(ptr)->set(r, c, *val);	
    return;
}

void AaronMatrix_integer_output_setRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void AaronMatrix_integer_output_setCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

void AaronMatrix_integer_output_setRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void AaronMatrix_integer_output_setCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronIntOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

// Indexed setters

void AaronMatrix_integer_output_setRowIndexed_integer(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::IntegerVector::iterator* out) {
    static_cast<AaronIntOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void AaronMatrix_integer_output_setColIndexed_integer(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::IntegerVector::iterator* out) {
    static_cast<AaronIntOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}

void AaronMatrix_integer_output_setRowIndexed_numeric(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::NumericVector::iterator* out) {
    static_cast<AaronIntOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void AaronMatrix_integer_output_setColIndexed_numeric(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::NumericVector::iterator* out) {
    static_cast<AaronIntOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}
