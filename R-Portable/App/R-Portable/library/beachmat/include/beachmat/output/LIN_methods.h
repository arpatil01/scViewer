#ifndef BEACHMAT_LIN_METHODS_WRITE_H
#define BEACHMAT_LIN_METHODS_WRITE_H

namespace beachmat { 

/****************************************
 * Defining the common output interface. 
 ****************************************/

// Getters:
template<typename T, class V>
void lin_output<T, V>::get_col(size_t c, Rcpp::IntegerVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::get_col(size_t c, Rcpp::NumericVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::get_row(size_t r, Rcpp::IntegerVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_output<T, V>::get_row(size_t r, Rcpp::NumericVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

// Setters:
template<typename T, class V>
void lin_output<T, V>::set_col(size_t c, Rcpp::IntegerVector::iterator out) {
    set_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::set_col(size_t c, Rcpp::NumericVector::iterator out) {
    set_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::set_row(size_t r, Rcpp::IntegerVector::iterator out) {
    set_row(r, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_output<T, V>::set_row(size_t r, Rcpp::NumericVector::iterator out) {
    set_row(r, out, 0, get_ncol());
    return;
}

/* Defining the general output interface. */ 

template<typename T, class V, class WTR>
general_lin_output<T, V, WTR>::general_lin_output(size_t nr, size_t nc) : writer(nr, nc) {}

// Getters:
template<typename T, class V, class WTR>
size_t general_lin_output<T, V, WTR>::get_nrow() const {
    return writer.get_nrow();
}

template<typename T, class V, class WTR>
size_t general_lin_output<T, V, WTR>::get_ncol() const {
    return writer.get_ncol();
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
T general_lin_output<T, V, WTR>::get(size_t r, size_t c) {
    return writer.get(r, c);
}

// Setters:
template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.set_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.set_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.set_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.set_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set(size_t r, size_t c, T in) {
    writer.set(r, c, in);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
    writer.set_col_indexed(c, N, idx, val);
    return; 
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
    writer.set_col_indexed(c, N, idx, val);
    return; 
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
    writer.set_row_indexed(r, N, idx, val);
    return; 
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
    writer.set_row_indexed(r, N, idx, val);
    return; 
}

// Other functions:
template<typename T, class V, class WTR>
Rcpp::RObject general_lin_output<T, V, WTR>::yield() {
    return writer.yield();
}

template<typename T, class V, class WTR>
std::unique_ptr<lin_output<T, V> > general_lin_output<T, V, WTR>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new general_lin_output<T, V, WTR>(*this));
}

}

#endif
