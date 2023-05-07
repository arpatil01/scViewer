#ifndef BEACHMAT_LIN_METHODS_READ_H
#define BEACHMAT_LIN_METHODS_READ_H

namespace beachmat { 

/****************************************
 * Defining the common input interface. 
 ****************************************/

// Basic getters.

template<typename T, class V>
void lin_matrix<T, V>::get_col(size_t c, Rcpp::IntegerVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_col(size_t c, Rcpp::NumericVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_row(size_t r, Rcpp::IntegerVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_row(size_t r, Rcpp::NumericVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

// Multi getters.

template<typename T, class V>
void lin_matrix<T, V>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out) {
    get_cols(it, n, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out) {
    get_cols(it, n, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out) {
    get_rows(it, n, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out) {
    get_rows(it, n, out, 0, get_ncol());
    return;
}

/* Defining the general interface. */

template<typename T, class V, class RDR>
general_lin_matrix<T, V, RDR>::general_lin_matrix(const Rcpp::RObject& incoming) : reader(incoming) {}

// Basic getters.

template<typename T, class V, class RDR>
size_t general_lin_matrix<T, V, RDR>::get_nrow() const {
    return reader.get_nrow();
}

template<typename T, class V, class RDR>
size_t general_lin_matrix<T, V, RDR>::get_ncol() const {
    return reader.get_ncol();
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class RDR>
T general_lin_matrix<T, V, RDR>::get(size_t r, size_t c) {
    return reader.get(r, c);
}

// Multi getters.

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_cols(it, n, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_cols(it, n, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_rows(it, n, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_rows(it, n, out, first, last);
    return;
}

// Other methods.

template<typename T, class V, class RDR>
std::unique_ptr<lin_matrix<T, V> > general_lin_matrix<T, V, RDR>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new general_lin_matrix<T, V, RDR>(*this));
}

template<typename T, class V, class RDR> 
Rcpp::RObject general_lin_matrix<T, V, RDR>::yield() const {
    return reader.yield();
}

}

#endif
