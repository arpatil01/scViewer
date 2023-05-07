#ifndef BEACHTEST_GET_VARSLICE_H
#define BEACHTEST_GET_VARSLICE_H
#include "beachtest.h"

template <class T, class M>  
Rcpp::List get_row_varslice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix cols) {
    if (cols.ncol()!=2) { 
        throw std::runtime_error("'cols' should be an integer matrix with two columns"); 
    }
    if (cols.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(cols)' should be equal to 'length(ordering)'");
    }
    Rcpp::List output(ordering.size());

    size_t r=0;
    for (auto o : ordering) {
        auto cur_bounds=cols.row(r);
        int left=cur_bounds[0]-1, right=cur_bounds[1];

        T target(right-left);
        ptr->get_row(o-1, target.begin(), left, right);
        output[r]=target;
        ++r;
    }

    return output;
}

template <class T, class M>  
Rcpp::List get_col_varslice (M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix rows) {
    if (rows.ncol()!=2) { 
        throw std::runtime_error("'rows' should be an integer matrix with two columns"); 
    }
    if (rows.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(rows)' should be equal to 'length(ordering)'");
    }
    Rcpp::List output(ordering.size());

    size_t c=0;
    for (auto o : ordering) {
        auto cur_bounds=rows.row(c);
        int left=cur_bounds[0]-1, right=cur_bounds[1];

        T target(right-left);
        ptr->get_col(o-1, target.begin(), left, right);
        output[c]=target;
        ++c;
    }

    return output;
}

#endif
