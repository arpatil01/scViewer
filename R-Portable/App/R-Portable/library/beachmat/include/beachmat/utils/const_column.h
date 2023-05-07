#ifndef BEACHMAT_CONST_COLUMN_H
#define BEACHMAT_CONST_COLUMN_H

#include "Rcpp.h"
#include "raw_structure.h"
#include <algorithm>

namespace beachmat {

/* A convenience class to obtain constant columns from a given matrix,
 * while supporting native sparse representations more-or-less transparently.
 */

template<class M>
class const_column {
public:
    const_column(M* mat, bool allow_sparse=true) : ref(mat), raws(mat->set_up_raw()), 
        Is_dense(mat->col_raw_type()=="dense"), Is_sparse(allow_sparse && mat->col_raw_type()=="sparse")
    {
        if (!Is_dense && !Is_sparse) {
            // repurposing the raw structure to hold some values.
            raws=raw_structure<typename M::vector>(mat->get_nrow()); 
        }
        return;
    }

    bool is_sparse () const { return Is_sparse; }

    bool is_dense () const { return Is_dense; }

    void fill(size_t c, size_t first, size_t last) {
        if (Is_dense || Is_sparse) {
            ref->get_col_raw(c, raws, first, last);
        } else {
            ref->get_col(c, raws.values.vec.begin(), first, last);
        }
        if (!Is_sparse) {
            raws.n=last - first;
            prev_start=first;
        }
        return;
    }

    void fill(size_t c) {
        fill(c, 0, ref->get_nrow());
        return;
    }

    size_t get_n () const {
        return raws.n;
    }

    // Not const, as Rcpp iterator conversions are problerefic.
    typename M::vector::iterator get_values() {
        if (!Is_dense && !Is_sparse) {
            return raws.values.vec.begin();
        } else {
            return raws.values_start;
        }
    }

    Rcpp::IntegerVector::iterator get_indices() {
        if (Is_sparse) {
            return raws.structure_start;
        }
        if (ref->get_nrow() > indices.size()) {
            indices=Rcpp::IntegerVector(ref->get_nrow());
            std::iota(indices.begin(), indices.end(), 0);
        }
        return indices.begin()+prev_start;
    }
private:
    M* ref;
    raw_structure<typename M::vector> raws;
    bool Is_dense, Is_sparse;

    Rcpp::IntegerVector indices; // deliberately copyable; values won't change, and reassignment won't matter.
    size_t prev_start=0;
};

}

#endif
