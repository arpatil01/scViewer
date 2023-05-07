#ifndef BEACHTEST_SET_ALL_H
#define BEACHTEST_SET_ALL_H
#include "beachtest.h"

template <class T, class O, class M>  // M is automatically deduced.
void set_row_all (O in, M ptr, Rcpp::IntegerVector ordering) {
    size_t r=0;
    T target(in->get_ncol());
    for (auto o : ordering) {
        in->get_row(o-1, target.begin());
        ptr->set_row(r, target.begin());
        ++r;
    }
    return;
}

template <class T, class O, class M>  // M is automatically deduced.
void set_col_all(O in, M ptr, Rcpp::IntegerVector ordering) {
    size_t c=0;
    T target(in->get_nrow());
    for (auto o : ordering) {
        in->get_col(o-1, target.begin());
        ptr->set_col(c, target.begin());
        ++c;
    }
    return;
}

template <class T, class O, class M>  // M is automatically deduced.
void set_single_all (O in, M ptr, Rcpp::IntegerVector rorder, Rcpp::IntegerVector corder) {
    size_t c=0;
    for (auto co : corder) {
        size_t r=0;
        for (auto ro : rorder) {
            ptr->set(r, c, in->get(ro-1, co-1));
            ++r;
        }
        ++c;
    }
    return;
}

#endif
