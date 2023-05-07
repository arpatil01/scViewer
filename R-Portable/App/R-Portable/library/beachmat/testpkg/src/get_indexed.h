#ifndef BEACHTEST_GET_INDEXED_H
#define BEACHTEST_GET_INDEXED_H

#include "beachtest.h"
#include "beachmat/utils/const_column.h"
#include <algorithm>

template <class O, class M>  // M is automatically deduced.
O get_indexed_all (M* ptr, Rcpp::IntegerVector ordering) {
    const size_t& nrows=ptr->get_nrow();
    O output(nrows, ordering.size());
    beachmat::const_column<M> col_holder(ptr);
    size_t c=0;

    for (auto o : ordering) {
        col_holder.fill(o-1);
        auto N=col_holder.get_n();
        auto idx=col_holder.get_indices();
        auto vals=col_holder.get_values();

        auto outcol=output.column(c);
        for (size_t i=0; i<N; ++i, ++idx, ++vals) {
            outcol[*idx] = *vals;
        }
        ++c;
    }
    return output;
}

template <class O, class M>  
O get_indexed_slice (M* ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int out_nrows=rend-rstart;
    O output(out_nrows, ordering.size());
    beachmat::const_column<M> col_holder(ptr);
    size_t c=0;

    for (auto o : ordering) {
        col_holder.fill(o-1, rstart, rend);
        auto N=col_holder.get_n();
        auto idx=col_holder.get_indices();
        auto vals=col_holder.get_values();

        auto curcol=output.column(c);
        for (size_t i=0; i<N; ++i, ++idx, ++vals) {
            curcol[*idx - rstart] = *vals;
        }
        ++c;
    }

    return output;
}

template <class M>  
Rcpp::List get_indexed_varslice (M* ptr, Rcpp::IntegerVector ordering, Rcpp::IntegerMatrix rows) {
    if (rows.ncol()!=2) { 
        throw std::runtime_error("'rows' should be an integer matrix with two columns"); 
    }
    if (rows.nrow()!=ordering.size()) {
        throw std::runtime_error("'nrow(rows)' should be equal to 'length(ordering)'");
    }
    Rcpp::List output(ordering.size());
    const size_t nrows=ptr->get_nrow();
    beachmat::const_column<M> col_holder(ptr);
    size_t c=0;

    for (auto o : ordering) {
        auto cur_bounds=rows.row(c);
        int left=cur_bounds[0]-1, right=cur_bounds[1];

        col_holder.fill(o-1, left, right);
        auto N=col_holder.get_n();
        auto idx=col_holder.get_indices();
        auto vals=col_holder.get_values();

        typename M::vector out(right-left);
        for (size_t i=0; i<N; ++i, ++idx, ++vals) {
            out[*idx - left] = *vals;
        }
        output[c]=out;
        ++c;
    }

    return output;
}

#endif
