#ifndef BEACHTEST_GET_ALL_H
#define BEACHTEST_GET_ALL_H
#include "beachtest.h"

template <class T, class O, class M>  // M is automatically deduced.
O get_row_all (M ptr, Rcpp::IntegerVector ordering) {
    const size_t& ncols=ptr->get_ncol();
    O output(ordering.size(), ncols);

	T target(ncols);
    size_t r=0;
    for (auto o : ordering){ 
        ptr->get_row(o-1, target.begin());
        auto outrow=output.row(r);
        std::copy(target.begin(), target.end(), outrow.begin());
        ++r;
    }
    return output;
}

template <class T, class O, class M>  // M is automatically deduced.
O get_col_all (M ptr, Rcpp::IntegerVector ordering) {
    const size_t& nrows=ptr->get_nrow();
    O output(nrows, ordering.size());

	T target(nrows);
    size_t c=0;
    for (auto o : ordering) {
        ptr->get_col(o-1, target.begin());
        auto outcol=output.column(c);
        std::copy(target.begin(), target.end(), outcol.begin());
        ++c;
    }
    return output;
}

template <class T, class O, class M>  // M is automatically deduced.
O get_single_all (M ptr, Rcpp::IntegerVector rorder, Rcpp::IntegerVector corder) {
    O output(rorder.size(), corder.size());
            
    size_t c=0;
    for (auto co : corder) {
        size_t i=c*rorder.size();
        for (auto ro : rorder) {
            output[i]=ptr->get(ro-1, co-1);
            ++i;
        }
        ++c;
    }

    return output;
}

#endif
