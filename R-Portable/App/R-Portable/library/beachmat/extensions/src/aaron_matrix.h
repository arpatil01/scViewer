#ifndef AARON_MATRIX_H
#define AARON_MATRIX_H

#include "Rcpp.h"

/* NOTE: input indices can be guaranteed to be within [0, nrow) for rows and [0, ncol) for columns.
 * For get_rows() and get_cols(), indices are also guaranteed to be strictly increasing.
 * Further checks are therefore unnecessary at this point.
 */

template<typename T, class V, class M>
class AaronMatrix {
public:
    AaronMatrix(const Rcpp::RObject& incoming) : mat(Rcpp::RObject(Rcpp::S4(incoming).slot("data"))), vec(mat) {}

    size_t get_nrow() const { return mat.nrow(); }
    size_t get_ncol() const { return mat.ncol(); }

    // Basic getters.
    T get(size_t r, size_t c) { return mat(r, c); }

    template<class Iter>
    void get_row(size_t r, Iter out, size_t first, size_t last) {
        auto currow=mat.row(r);
        std::copy(currow.begin()+first, currow.begin()+last, out);
        return;
    }

    template<class Iter>
    void get_col(size_t c, Iter out, size_t first, size_t last) {
        auto curcol=mat.column(c);
        std::copy(curcol.begin()+first, curcol.begin()+last, out);
        return;
    }
    
    // Multi getters.
    template<class Iter>
    void get_rows(Rcpp::IntegerVector::iterator r, size_t n, Iter out, size_t first, size_t last) {
        for (size_t c=first; c<last; ++c) {
            auto curcol=mat.column(c);
            auto it=curcol.begin();
            auto r_copy=r;
            for (size_t i=0; i<n; ++i, ++out, ++r_copy) {  
                (*out)=*(it + *r_copy);
            }
        }
    }

    template<class Iter>
    void get_cols(Rcpp::IntegerVector::iterator c, size_t n, Iter out, size_t first, size_t last) {
        for (size_t i=0; i<n; ++i) {
            get_col(*c, out, first, last);
            out+=last - first;
            ++c;
        }
        return;
    }

private:
    M mat;
    V vec; // costless conversion to a vector.
    Rcpp::IntegerVector indices;
};

#endif
