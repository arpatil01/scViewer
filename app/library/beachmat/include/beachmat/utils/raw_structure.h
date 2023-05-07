#ifndef BEACHMAT_RAW_STRUCTURE_H
#define BEACHMAT_RAW_STRUCTURE_H

#include "Rcpp.h"
#include "copyable_vector.h"

namespace beachmat {

/* This class is strictly a convenience class to hold various data members together:
 *
 * - values_start may or may not point to any valid array.
 * - structure_start may or may not point to any valid array.
 * - values may or may not contain an empty vector. 
 * - structure may or may not contain an empty vector. 
 *
 * What is valid or not depends on the specific matrix representation.
 */

template<class V>
struct raw_structure {
    raw_structure(size_t nv=0, size_t ns=0) : values(nv), structure(ns) {}

    size_t n=0;
    copyable_holder<V> values;
    typename V::iterator values_start;
    copyable_holder<Rcpp::IntegerVector> structure;
    typename Rcpp::IntegerVector::iterator structure_start;
};

}

#endif
