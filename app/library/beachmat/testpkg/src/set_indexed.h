#ifndef BEACHTEST_SET_INDEXED_H
#define BEACHTEST_SET_INDEXED_H
#include "beachtest.h"

template <class V, class OX>  
Rcpp::RObject set_row_indexed (OX optr, const Rcpp::IntegerVector& indices, Rcpp::List subindices) {
    if (indices.size()!=subindices.size()) { 
        throw std::runtime_error("'index' and 'subindices' should be of the same length");
    }

    for (size_t i=0; i<indices.size(); ++i) {
        Rcpp::List cursub=subindices[i];
        if (cursub.size()!=2) { 
            throw std::runtime_error("each element of 'subindices' should be a list of length 2");
        }

        Rcpp::IntegerVector I=cursub[0];
        V X=cursub[1];
        if (X.size()!=I.size()) {
            throw std::runtime_error("length of vectors in each element of 'subindices' should be equal");
        }
        
        // Some work required to get to zero-based indexing.
        Rcpp::IntegerVector I2(Rcpp::clone(I));
        for (auto& x : I2) { --x; } 

        optr->set_row_indexed(indices[i] - 1, I2.size(), I2.begin(), X.begin()); 
    }
    return optr->yield();
}

template <class V, class OX>  
Rcpp::RObject set_col_indexed (OX optr, const Rcpp::IntegerVector& indices, Rcpp::List subindices) {
    if (indices.size()!=subindices.size()) { 
        throw std::runtime_error("'index' and 'subindices' should be of the same length");
    }

    for (size_t i=0; i<indices.size(); ++i) {
        Rcpp::List cursub=subindices[i];
        if (cursub.size()!=2) { 
            throw std::runtime_error("each element of 'subindices' should be a list of length 2");
        }

        Rcpp::IntegerVector I=cursub[0];
        V X=cursub[1];
        if (X.size()!=I.size()) {
            throw std::runtime_error("length of vectors in each element of 'subindices' should be equal");
        }
        
        // Some work required to get to zero-based indexing.
        Rcpp::IntegerVector I2(Rcpp::clone(I));
        for (auto& x : I2) { --x; } 

        optr->set_col_indexed(indices[i] - 1, I2.size(), I2.begin(), X.begin()); 
    }

    return optr->yield();
}

#endif
