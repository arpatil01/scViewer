#ifndef BEACHTEST_GET_SLICE_H
#define BEACHTEST_GET_SLICE_H
#include "beachtest.h"

template <class T, class O, class M>  
O get_row_slice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector cols) {
    if (cols.size()!=2) { 
        throw std::runtime_error("'cols' should be an integer vector of length 2"); 
    }
    const int cstart=cols[0]-1, cend=cols[1];
    const int ncols=cend-cstart;    

    O output(ordering.size(), ncols);
    T target(ncols);
    size_t r=0;
    for (auto o : ordering) {
        ptr->get_row(o-1, target.begin(), cstart, cend);
        auto currow=output.row(r);
        std::copy(target.begin(), target.end(), currow.begin());
        ++r;
    }

    return output;
}

template <class T, class O, class M>  
O get_col_slice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int nrows=rend-rstart;    

    O output(nrows, ordering.size());
    T target(nrows);
    size_t c=0;
    for (auto o : ordering) {
        ptr->get_col(o-1, target.begin(), rstart, rend);
        auto curcol=output.column(c);
        std::copy(target.begin(), target.end(), curcol.begin());
        ++c;
    }

    return output;
}

#endif
