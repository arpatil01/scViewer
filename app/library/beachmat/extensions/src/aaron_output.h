#ifndef AARON_MATRIX_H
#define AARON_MATRIX_H

#include "Rcpp.h"

template<typename T, class V, class M>
class AaronOutput {
public:
    AaronOutput(size_t nr, size_t nc) : mat(nr, nc), vec(mat) {}

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
    
    // Basic setters.
    void set(size_t r, size_t c, T val) { 
        *(mat.begin() + r + c*mat.nrow())=val; 
        return;
    }

    template<class Iter>
    void set_row(size_t r, Iter out, size_t first, size_t last) {
        auto currow=mat.row(r);
        std::copy(out, out+last-first, currow.begin()+first);
        return;
    }

    template<class Iter>
    void set_col(size_t c, Iter out, size_t first, size_t last) {
        auto curcol=mat.column(c);
        std::copy(out, out+last-first, curcol.begin()+first);
        return;
    }
 
    // Indexed setters.
    template<class Iter>
    void set_row_indexed(size_t r, size_t n, Rcpp::IntegerVector::iterator idx, Iter out) {
        auto currow=mat.row(r);
        for (size_t i=0; i<n; ++i, ++idx, ++out) { currow[*idx]=*out; }
        return;
    }

    template<class Iter>
    void set_col_indexed(size_t c, size_t n, Rcpp::IntegerVector::iterator idx, Iter out) {
        auto curcol=mat.column(c);
        for (size_t i=0; i<n; ++i, ++idx, ++out) { curcol[*idx]=*out; }
        return;
    }
  
    Rcpp::RObject yield () const {
        Rcpp::S4 out("AaronMatrix");
        out.slot("data")=mat;
        return Rcpp::RObject(out);
    }
private:
    M mat;
    V vec; // costless conversion to a vector.
};

#endif
