#ifndef BEACHMAT_SIMPLE_READER_H
#define BEACHMAT_SIMPLE_READER_H

#include "Rcpp.h"

#include "../utils/utils.h"
#include "../utils/dim_checker.h"
#include "../utils/raw_structure.h"

#include <stdexcept>
#include <algorithm>

namespace beachmat {

/*** Class definition ***/

template<typename T, class V>
class simple_reader : public dim_checker {
public:    
    simple_reader(const Rcpp::RObject&);
    ~simple_reader() = default;
    simple_reader(const simple_reader&) = default;
    simple_reader& operator=(const simple_reader&) = default;
    simple_reader(simple_reader&&) = default;
    simple_reader& operator=(simple_reader&&) = default;

    // Basic getters.
    T get(size_t, size_t);

    template <class Iter>
    void get_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    // Specialized getters.
    raw_structure<V> set_up_raw () const {
        return raw_structure<V>();
    }

    void get_col_raw(size_t c, raw_structure<V>& in, size_t first, size_t last) {
        check_colargs(c, first, last);
        in.values_start=get_const_col(c, first);
        return;
    }

    void get_row_raw(size_t r, raw_structure<V>& in, size_t first, size_t last) {
        check_rowargs(r, first, last);
        return;
    }

    static std::string col_raw_type () { return "dense"; }

    static std::string row_raw_type () { return "none"; }

    // Multi getters.
    template <class Iter>
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    // Miscellaneous.
    Rcpp::RObject yield() const { return original; }

    static std::string get_class() { return "matrix"; }

    static std::string get_package() { return "base"; }
private:
    Rcpp::RObject original;
    V mat;

    typename V::iterator get_const_col(size_t c, size_t first) {
        return mat.begin() + first + c*(this->nrow);
    }
};

/*** Constructor definitions ***/

template<typename T, class V>
simple_reader<T, V>::simple_reader(const Rcpp::RObject& incoming) : original(incoming) { 
    if (!incoming.hasAttribute("dim")) { 
        throw std::runtime_error("matrix object should have 'dim' attribute"); 
    }
    this->fill_dims(incoming.attr("dim"));
    const size_t& NC=this->ncol; 

    if (incoming.sexp_type()!=mat.sexp_type()) { 
        throw std::runtime_error(std::string("matrix should be ") + translate_type(mat.sexp_type()));
    }
    mat=incoming;
    if (static_cast<size_t>(mat.size())!=(this->nrow)*NC) {
        throw std::runtime_error("length of matrix is inconsistent with its dimensions"); 
    }
    return;
}

/*** Basic getter methods ***/

template<typename T, class V>
T simple_reader<T, V>::get(size_t r, size_t c) { 
    check_oneargs(r, c);
    return mat[c*(this->nrow)+r]; // do NOT use get_const_col(c, r), as it fails for Strings.
}

template<typename T, class V>
template<class Iter>
void simple_reader<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    const size_t& NR=this->nrow;
    auto src=get_const_col(first, r);
    for (size_t col=first; col<last; ++col, src+=NR, ++out) { (*out)=(*src); }
    return;
}

template<typename T, class V>
template<class Iter>
void simple_reader<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    auto src=get_const_col(c, 0);
    std::copy(src+first, src+last, out);
    return;
}

/*** Multi getter methods ***/

template<typename T, class V>
template<class Iter>
void simple_reader<T, V>::get_rows(Rcpp::IntegerVector::iterator rIt, size_t n, Iter out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);

    for (size_t c=first; c<last; ++c) {
        auto it=get_const_col(c, 0);
        auto rIt_copy=rIt;
        for (size_t i=0; i<n; ++i, ++out, ++rIt_copy) {  
            (*out)=*(it + *rIt_copy);
        }
    }
    return;
}

template<typename T, class V>
template<class Iter>
void simple_reader<T, V>::get_cols(Rcpp::IntegerVector::iterator cIt, size_t n, Iter out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);
    size_t nrows=last - first;
    for (size_t i=0; i<n; ++i, ++cIt, out+=nrows) {
        get_col(*cIt, out, first, last);
    }
    return;
}

}

#endif
