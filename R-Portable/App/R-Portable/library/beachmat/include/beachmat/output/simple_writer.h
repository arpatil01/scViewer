#ifndef BEACHMAT_SIMPLE_WRITER_H
#define BEACHMAT_SIMPLE_WRITER_H

#include "Rcpp.h"

#include "../utils/utils.h"
#include "../utils/dim_checker.h"

#include <algorithm>

namespace beachmat {

/*** Class definition ***/

template<typename T, class V>
class simple_writer : public dim_checker {
public:
    simple_writer(size_t, size_t);
    ~simple_writer() = default;
    simple_writer(const simple_writer&) = default;
    simple_writer& operator=(const simple_writer&) = default;
    simple_writer(simple_writer&&) = default;
    simple_writer& operator=(simple_writer&&) = default;

    // Setters:
    template <class Iter>
    void set_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void set_col(size_t, Iter, size_t, size_t);

    void set(size_t, size_t, T);

    template <class Iter>
    void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Iter);

    template <class Iter>
    void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Iter);

    // Getters:
    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_row(size_t, Iter, size_t, size_t);

    T get(size_t, size_t);
    
    // Other:
    Rcpp::RObject yield();

    static std::string get_class() { return "matrix"; }

    static std::string get_package() { return "base"; }
private:
    V data;    
};

/*** Constructor definition ***/

template<typename T, class V>
simple_writer<T, V>::simple_writer(size_t nr, size_t nc) : dim_checker(nr, nc) { 
    (this->data)=V(nr*nc);
    return; 
}

/*** Setter methods ***/

template<typename T, class V>
template<class Iter>
void simple_writer<T, V>::set_col(size_t c, Iter in, size_t start, size_t end) {
    check_colargs(c, start, end);
    std::copy(in, in + end - start, data.begin()+c*(this->nrow)+start); 
    return;
}

template<typename T, class V>
template<class Iter>
void simple_writer<T, V>::set_row(size_t r, Iter in, size_t start, size_t end) {
    check_rowargs(r, start, end);
    const size_t& NR=this->nrow;
    auto mIt=data.begin() + r + start * NR;
    for (size_t c=start; c<end; ++c, mIt+=NR, ++in) {
        (*mIt)=*in;        
    }
    return;
}

template<typename T, class V>
void simple_writer<T, V>::set(size_t r, size_t c, T in) {
    check_oneargs(r, c);
    data[r + (this->nrow)*c]=in;
    return;
}

template<typename T, class V>
template <class Iter>
void simple_writer<T, V>::set_col_indexed(size_t c, size_t n, Rcpp::IntegerVector::iterator idx, Iter in) {
    check_colargs(c);
    auto current=data.begin() + c * (this->nrow);
    for (size_t i=0; i<n; ++i, ++idx, ++in) {
        *(current + *idx) = *in;
    }
    return;
}

template<typename T, class V>
template <class Iter>
void simple_writer<T, V>::set_row_indexed(size_t r, size_t n, Rcpp::IntegerVector::iterator idx, Iter in) {
    check_rowargs(r);
    auto current=data.begin() + r;
    for (size_t i=0; i<n; ++i, ++idx, ++in) {
        *(current + (*idx)*(this->nrow)) = *in;
    }
    return;
}

/*** Getter methods ***/

template<typename T, class V>
template<class Iter>
void simple_writer<T, V>::get_row(size_t r, Iter out, size_t start, size_t end) {
    check_rowargs(r, start, end);
    const size_t& NR=this->nrow;
    auto src=data.begin()+start*NR+r;
    for (size_t col=start; col<end; ++col, src+=NR, ++out) { (*out)=(*src); }
    return;
}

template<typename T, class V>
template<class Iter>
void simple_writer<T, V>::get_col(size_t c, Iter out, size_t start, size_t end) {
    check_colargs(c, start, end);
    auto src=data.begin() + c*(this->nrow);
    std::copy(src+start, src+end, out);
    return;
}

template<typename T, class V>
T simple_writer<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    return data[c*(this->nrow)+r];
}

/*** Output function ***/

template<typename T, class V>
Rcpp::RObject simple_writer<T, V>::yield() {
    Rcpp::RObject out(SEXP(this->data));
    out.attr("dim") = Rcpp::IntegerVector::create(this->nrow, this->ncol); 
    return out;
}

}

#endif
