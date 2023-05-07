#ifndef BEACHTEST_SET_SLICE_H
#define BEACHTEST_SET_SLICE_H
#include "beachtest.h"

template <class T, class O, class M>  // M is automatically deduced.
void set_row_slice (O in, M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector cols) {
    if (cols.size()!=2) { 
        throw std::runtime_error("'cols' should be an integer vector of length 2"); 
    }
    const int cstart=cols[0]-1, cend=cols[1];
    const int ncols=cend-cstart;    

    size_t r=0;
    T target(ncols);
    for (auto o : ordering) {
        in->get_row(o-1, target.begin(), cstart, cend);
        ptr->set_row(r, target.begin(), cstart, cend);
        ++r;
    }
    return;
}

template <class T, class O, class M>  // M is automatically deduced.
void set_col_slice(O in, M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int nrows=rend-rstart;

    size_t c=0;
    T target(nrows);
    for (auto o : ordering) {
        in->get_col(o-1, target.begin(), rstart, rend);
        ptr->set_col(c, target.begin(), rstart, rend);
        ++c;
    }
    return;
}

#endif
