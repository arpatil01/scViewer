#ifndef BEACHTEST_SET_VARSLICE_H
#define BEACHTEST_SET_VARSLICE_H
#include "beachtest.h"

template <class T, class O, class M>  
void set_row_varslice (O in, M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix cols) {
    if (cols.ncol()!=2) { 
        throw std::runtime_error("'cols' should be an integer matrix with two columns"); 
    }
    if (cols.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(cols)' should be equal to 'length(ordering)'");
    }

    size_t r=0;
    for (auto o : ordering) {
        auto cur_bounds=cols.row(r);
        int left=cur_bounds[0]-1, right=cur_bounds[1];

        T target(right-left);
        in->get_row(o-1, target.begin(), left, right);
        ptr->set_row(r, target.begin(), left, right);
        ++r;
    }

    return;
}

template <class T, class O, class M>  
void set_col_varslice (O in, M ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix rows) {
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
        in->get_col(o-1, target.begin(), left, right);
        ptr->set_col(c, target.begin(), left, right);
        ++c;
    }

    return;
}

#endif
