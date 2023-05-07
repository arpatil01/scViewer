#include "aaron_output.h"
#include "exports.h"

// Demonstrating with integer matrices.

typedef AaronOutput<Rcpp::String, Rcpp::StringVector, Rcpp::StringMatrix> AaronStrOut;

// Constructor, destructors and clones.

void * AaronMatrix_character_output_create(size_t nr, size_t nc) {
    return static_cast<void*>(new AaronStrOut(nr, nc));
}

void AaronMatrix_character_output_destroy(void * ptr) {
    delete static_cast<AaronStrOut*>(ptr);
    return;
}

void * AaronMatrix_character_output_clone(void * ptr) {
    AaronStrOut* old=static_cast<AaronStrOut*>(ptr);
    return static_cast<void*>(new AaronStrOut(*old));
}

SEXP AaronMatrix_character_output_yield(void * ptr) {
    return static_cast<AaronStrOut*>(ptr)->yield();
}

// Basic getters

void AaronMatrix_character_output_get(void * ptr, size_t r, size_t c, Rcpp::String* val) {
    *val=static_cast<AaronStrOut*>(ptr)->get(r, c);	
    return;
}

void AaronMatrix_character_output_getRow(void * ptr, size_t r, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void AaronMatrix_character_output_getCol(void * ptr, size_t c, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Basic setters

void AaronMatrix_character_output_set(void * ptr, size_t r, size_t c, Rcpp::String* val) {
    static_cast<AaronStrOut*>(ptr)->set(r, c, *val);	
    return;
}

void AaronMatrix_character_output_setRow(void * ptr, size_t r, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void AaronMatrix_character_output_setCol(void * ptr, size_t c, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

// Indexed setters

void AaronMatrix_character_output_setRowIndexed(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::StringVector::iterator* out) {
    static_cast<AaronStrOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void AaronMatrix_character_output_setColIndexed(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::StringVector::iterator* out) {
    static_cast<AaronStrOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}

