#ifndef BEACHTEST_GET_MULTI_H
#define BEACHTEST_GET_MULTI_H
#include "beachtest.h"

template <class T, class O, class M>  
O get_multirow_all (M ptr, Rcpp::IntegerVector chosen) {
    O output(chosen.size(), ptr->get_ncol());
    T target(output);

    // For sanity's sake, all indices are zero-indexed in .Call().
    ptr->get_rows(chosen.begin(), chosen.size(), target.begin());
    return output;
}

template <class T, class O, class M>  
O get_multirow_slice (M ptr, Rcpp::IntegerVector chosen, Rcpp::IntegerVector cols) {
    if (cols.size()!=2) { 
        throw std::runtime_error("'cols' should be an integer vector of length 2"); 
    }
    const int cstart=cols[0]-1, cend=cols[1];
    const int ncols=cend-cstart;    

    O output(chosen.size(), ncols);
    T target(output);
    ptr->get_rows(chosen.begin(), chosen.size(), target.begin(), cstart, cend);
    return output;
}

template <class T, class O, class M>  
O get_multicol_all (M ptr, Rcpp::IntegerVector chosen) {
    O output(ptr->get_nrow(), chosen.size());
    T target(output);
    ptr->get_cols(chosen.begin(), chosen.size(), target.begin());
    return output;
}

template <class T, class O, class M>  
O get_multicol_slice (M ptr, Rcpp::IntegerVector chosen, Rcpp::IntegerVector rows) {
    if (rows.size()!=2) { 
        throw std::runtime_error("'rows' should be an integer vector of length 2"); 
    }
    const int rstart=rows[0]-1, rend=rows[1];
    const int nrows=rend-rstart;    

    O output(nrows, chosen.size());
    T target(output);
    ptr->get_cols(chosen.begin(), chosen.size(), target.begin(), rstart, rend);
    return output;
}

#endif
