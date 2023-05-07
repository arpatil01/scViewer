#include "aaron_matrix.h"
#include "exports.h"

// Demonstrating with integer matrices.

typedef AaronMatrix<Rcpp::String, Rcpp::StringVector, Rcpp::StringMatrix> AaronStrMat;

// Constructor, destructors and clones.

void * AaronMatrix_character_input_create(SEXP incoming) {
    return static_cast<void*>(new AaronStrMat(incoming));
}

void AaronMatrix_character_input_destroy(void * ptr) {
    delete static_cast<AaronStrMat*>(ptr);
    return;
}

void * AaronMatrix_character_input_clone(void * ptr) {
    AaronStrMat* old=static_cast<AaronStrMat*>(ptr);
    return static_cast<void*>(new AaronStrMat(*old));
}

// Basic getters

void AaronMatrix_character_input_dim(void* ptr, size_t* nr, size_t* nc){ 
    AaronStrMat* thing=static_cast<AaronStrMat*>(ptr);
    *nr=thing->get_nrow();
    *nc=thing->get_ncol();
    return;
}

void AaronMatrix_character_input_get(void * ptr, size_t r, size_t c, Rcpp::String* val) {
    *val=static_cast<AaronStrMat*>(ptr)->get(r, c);	
    return;
}

void AaronMatrix_character_input_getRow(void * ptr, size_t r, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrMat*>(ptr)->get_row(r, *out, first, last);
    return;
}

void AaronMatrix_character_input_getCol(void * ptr, size_t c, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrMat*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Multi getters

void AaronMatrix_character_input_getRows(void * ptr, Rcpp::IntegerVector::iterator* r, size_t n, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrMat*>(ptr)->get_rows(*r, n, *out, first, last);
    return;
}

void AaronMatrix_character_input_getCols(void * ptr, Rcpp::IntegerVector::iterator* c, size_t n, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<AaronStrMat*>(ptr)->get_cols(*c, n, *out, first, last);
    return;
}

