#ifndef EXPORTS_H
#define EXPORTS_H
#include "Rcpp.h"

extern "C" {

void * AaronMatrix_character_input_create(SEXP);

void AaronMatrix_character_input_destroy(void *);

void * AaronMatrix_character_input_clone(void *);

void AaronMatrix_character_input_dim(void*, size_t*, size_t*);
 
void AaronMatrix_character_input_get(void *, size_t, size_t, Rcpp::String*);

void AaronMatrix_character_input_getRow(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getCol(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getRows(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getCols(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void * AaronMatrix_character_output_create(size_t, size_t);

void AaronMatrix_character_output_destroy(void *);

void * AaronMatrix_character_output_clone(void *);

SEXP AaronMatrix_character_output_yield(void *);

void AaronMatrix_character_output_get(void *, size_t, size_t, Rcpp::String*);

void AaronMatrix_character_output_getRow(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_output_getCol(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_output_set(void *, size_t, size_t, Rcpp::String*);

void AaronMatrix_character_output_setRow(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_output_setCol(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_output_setRowIndexed(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::StringVector::iterator*);

void AaronMatrix_character_output_setColIndexed(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::StringVector::iterator*);

void * AaronMatrix_integer_input_create (SEXP);

void AaronMatrix_integer_input_destroy (void *);

void * AaronMatrix_integer_input_clone (void *);

void AaronMatrix_integer_input_dim(void*, size_t*, size_t*);
 
void AaronMatrix_integer_input_get(void *, size_t, size_t, int*);

void AaronMatrix_integer_input_getRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getRows_integer(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCols_integer(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getRows_numeric(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCols_numeric(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void * AaronMatrix_integer_output_create (size_t, size_t);

void AaronMatrix_integer_output_destroy (void *);

void * AaronMatrix_integer_output_clone (void *);

SEXP AaronMatrix_integer_output_yield (void *);

void AaronMatrix_integer_output_get(void *, size_t, size_t, int*);

void AaronMatrix_integer_output_getRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_getCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_getRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_getCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_set(void *, size_t, size_t, int*);

void AaronMatrix_integer_output_setRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_setCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_setRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_setCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_output_setRowIndexed_integer(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::IntegerVector::iterator*);

void AaronMatrix_integer_output_setColIndexed_integer(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::IntegerVector::iterator*);

void AaronMatrix_integer_output_setRowIndexed_numeric(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::NumericVector::iterator*);

void AaronMatrix_integer_output_setColIndexed_numeric(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::NumericVector::iterator*);

}

#endif
