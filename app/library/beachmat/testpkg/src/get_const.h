#ifndef BEACHTEST_GET_CONST_H
#define BEACHTEST_GET_CONST_H

#include "beachtest.h"
#include "beachmat/utils/const_column.h"

template <class O, class M>  // M is automatically deduced.
O get_const_all (M* ptr, Rcpp::IntegerVector ordering) {
    const size_t nrows=ptr->get_nrow();
    O output(nrows, ordering.size());
    beachmat::const_column<M> col_holder(ptr, false);
    size_t c=0;

    for (auto o : ordering) {
        col_holder.fill(o-1);
        auto val=col_holder.get_values();
        auto outcol=output.column(c);
        std::copy(val, val + nrows, outcol.begin());
        ++c;
    }
    return output;
}

template <class O, class M>  
O get_const_slice (M* ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int nrows=rend-rstart;    
    O output(nrows, ordering.size());
    beachmat::const_column<M> col_holder(ptr, false);
    size_t c=0;

    for (auto o : ordering) {
        col_holder.fill(o-1, rstart, rend);
        auto val=col_holder.get_values();
        auto outcol=output.column(c);
        std::copy(val, val + nrows, outcol.begin());
        ++c;
    }

    return output;
}

template <class M>  
Rcpp::List get_const_varslice (M* ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix rows) {
    if (rows.ncol()!=2) { 
        throw std::runtime_error("'rows' should be an integer matrix with two columns"); 
    }
    if (rows.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(rows)' should be equal to 'length(ordering)'");
    }
    Rcpp::List output(ordering.size());
    beachmat::const_column<M> col_holder(ptr, false);
    size_t c=0;

    for (auto o : ordering) {
        auto cur_bounds=rows.row(c);
        int left=cur_bounds[0]-1, right=cur_bounds[1];
        col_holder.fill(o-1, left, right);
        auto val=col_holder.get_values();

        typename M::vector out(right-left);
        std::copy(val, val + out.size(), out.begin());
        output[c]=out;
        ++c;
    }

    return output;
}

#endif
