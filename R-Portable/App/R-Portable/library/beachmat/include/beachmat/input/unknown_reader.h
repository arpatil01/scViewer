#ifndef UNKNOWN_READER_H
#define UNKNOWN_READER_H

#include "Rcpp.h"

#include "../utils/dim_checker.h"
#include "../utils/copyable_vector.h"
#include "../utils/raw_structure.h"

#include <algorithm>

namespace beachmat {

/* The 'unknown_reader' class will realize chunks of the input RObject
 * upon request from any calling function. This was designed for 
 * DelayedMatrix objects, to avoid reimplementing arbitrary delayed 
 * R operations in C++; however, it is also useful for unknown matrices.
 */

template<typename T, class V>
class unknown_reader : public dim_checker {
public:    
    unknown_reader(const Rcpp::RObject&);
    ~unknown_reader() = default;
    unknown_reader(const unknown_reader&) = default;
    unknown_reader& operator=(const unknown_reader&) = default;
    unknown_reader(unknown_reader&&) = default;
    unknown_reader& operator=(unknown_reader&&) = default;

    T get(size_t, size_t);

    template <class Iter> 
    void get_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    // Specialized getters.
    raw_structure<V> set_up_raw () const {
        return raw_structure<V>();
    }

    void get_col_raw(size_t c, raw_structure<V>& in, size_t first, size_t last) {
        check_colargs(c, first, last);
        return;
    }

    void get_row_raw(size_t r, raw_structure<V>& in, size_t first, size_t last) {
        check_rowargs(r, first, last);
        return;
    }

    static std::string col_raw_type () { return "none"; }

    static std::string row_raw_type () { return "none"; }

    // Multi getters.
    template <class Iter>
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    // Miscellaneous.
    Rcpp::RObject yield() const { return original; }

    std::string get_class() const { return ""; }

    std::string get_package() const { return ""; }
private:
    Rcpp::RObject original;
    Rcpp::Environment beachenv;
    Rcpp::Function realizer;
    
    V storage; // this does not need to be copyable, as the entire object is replaced when a new block is loaded.
    void update_storage_by_row(size_t, size_t, size_t);
    void update_storage_by_col(size_t, size_t, size_t);
    
    size_t storage_start_row, storage_end_row, storage_start_col, storage_end_col; 
    bool oncol;

    Rcpp::IntegerVector row_chunk_bounds, col_chunk_bounds; // Read-only, does not need to be copyable.
    size_t chunk_index;

    static bool reload_chunk (size_t primary_value, size_t& primary_start, size_t& primary_end, 
        size_t& primary_bound_index, const Rcpp::IntegerVector& primary_bounds, 
        size_t secondary_val_1, size_t secondary_val_2, size_t& secondary_start, size_t& secondary_end); 

    // Values to be passed to R functions for realization - entries are modifiable and so this needs to be copiable.
    copyable_holder<Rcpp::IntegerVector> indices, slices;
    copyable_holder<Rcpp::LogicalVector> do_transpose;
};

/* Constructor definition. */

template<typename T, class V>
unknown_reader<T, V>::unknown_reader(const Rcpp::RObject& in) : original(in), 
        beachenv(Rcpp::Environment::namespace_env("beachmat")), realizer(beachenv["realizeByRange"]), 
        storage_start_row(0), storage_end_row(0), storage_start_col(0), storage_end_col(0), 
        oncol(false), 
        chunk_index(0),
        indices(2), slices(2), do_transpose(1) { 

    Rcpp::Function getdims(beachenv["setupUnknownMatrix"]);
    Rcpp::List dimdata=getdims(in);

    Rcpp::IntegerVector matdims(dimdata[0]);
    this->fill_dims(matdims);

    row_chunk_bounds=Rcpp::IntegerVector(dimdata[1]);
    col_chunk_bounds=Rcpp::IntegerVector(dimdata[2]);

    do_transpose.vec[0]=1;
    return;
}

/* Define storage-related methods. */

template<typename T, class V>
bool unknown_reader<T, V>::reload_chunk (size_t primary_value, size_t& primary_start, size_t& primary_end, 
        size_t& primary_bound_index, const Rcpp::IntegerVector& primary_bounds, 
        size_t secondary_val_1, size_t secondary_val_2, size_t& secondary_start, size_t& secondary_end) { 

    // We assume that all inputs are valid, as check_* functions are called upstream.
    const bool is_below = (primary_value < primary_start),
        is_above = (primary_value >= primary_end);
    const bool reobtain_primary=(is_below || is_above);

    // Determining which chunks to obtain on the primary dimension.
    if (reobtain_primary) {
        const int pval=primary_value;

        if (is_below) { 
            /* primary_bound_index always points at the end boundary of the current chunk.
             * In this scope, primary_bound_index CANNOT be less than 2, 
             * as is_below cannot be true if we're at the first chunk.
             */
            --primary_bound_index;
            if (primary_bounds[primary_bound_index-1] > pval) { 
                primary_bound_index=std::upper_bound(primary_bounds.begin() + 1, primary_bounds.begin() + primary_bound_index, pval) - primary_bounds.begin();
            }
            
        } else { 
            /* Note that the first search always ends up here as 'primary_start=primary_end=0'.
             * This means that 'is_below' must be false while 'is_above' must be true. 
             * We then increment primary_bound_index from 0 to 1, getting us to the first chunk.
             */
            ++primary_bound_index;
            if (primary_bounds[primary_bound_index] <= pval) {
                primary_bound_index=std::upper_bound(primary_bounds.begin() + primary_bound_index + 1, primary_bounds.end(), pval) - primary_bounds.begin();
            }
        }

        primary_end=primary_bounds[primary_bound_index];
        primary_start=primary_bounds[primary_bound_index-1];

    } else if (secondary_val_1 >= secondary_start && secondary_val_2 <= secondary_end) {
        // If the requested values are a subset of those already cached, there's no need to reload the chunk.
        return false;
    }

    secondary_start=secondary_val_1;
    secondary_end=secondary_val_2;
    return true;
}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_row(size_t r, size_t first, size_t last) {
    if (oncol) { // reset values to force 'reload_chunk' to give 'true'.
        storage_start_row=0;
        storage_end_row=0;
        chunk_index=0;
        oncol=false;
    }
        
    if (reload_chunk(r, storage_start_row, storage_end_row,
                chunk_index, row_chunk_bounds,
                first, last, storage_start_col, storage_end_col)) {

        indices.vec[0] = storage_start_row; 
        indices.vec[1] = storage_end_row - storage_start_row;
        slices.vec[0] = storage_start_col;
        slices.vec[1] = storage_end_col - storage_start_col;
        storage = realizer(original, indices.vec, slices.vec, do_transpose.vec); // Transposed, so storage is effectively row-major!
    }
    return;
}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_col(size_t c, size_t first, size_t last) {
    if (!oncol) { // reset values to force 'reload_chunk' to give 'true'.
        storage_start_col=0;
        storage_end_col=0;
        chunk_index=0;
        oncol=true;
    }

    if (reload_chunk(c, storage_start_col, storage_end_col,
                chunk_index, col_chunk_bounds, 
                first, last, storage_start_row, storage_end_row)) {

        indices.vec[0] = storage_start_col;
        indices.vec[1] = storage_end_col - storage_start_col;
        slices.vec[0] = storage_start_row; 
        slices.vec[1] = storage_end_row - storage_start_row;
        storage = realizer(original, slices.vec, indices.vec);
    }
    return;
}

/*** Basic getter methods ***/

template<typename T, class V>
T unknown_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    update_storage_by_col(c, 0, this->nrow); // keeps the whole block for further 'get()' queries.
    return storage[(c - storage_start_col) * this->nrow + r];
}

template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    update_storage_by_row(r, first, last);

    // It's effectively row-major storage due the transposition.
    auto src=storage.begin() + 
        (r - storage_start_row) * (storage_end_col - storage_start_col) +
        (first - storage_start_col);
    std::copy(src, src + (last - first), out);
    return;
}
 
template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    update_storage_by_col(c, first, last);

    auto src=storage.begin() + 
        (c - storage_start_col) * (storage_end_row - storage_start_row) + 
        (first - storage_start_row);
    std::copy(src, src + (last - first), out);
    return;
}

/*** Multi getter methods ***/

template<typename T, class V>
template<class Iter>
void unknown_reader<T, V>::get_rows(Rcpp::IntegerVector::iterator rIt, size_t n, Iter out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);

    // Need to make a copy of the indexed (1-indexed) to pass to the function.
    Rcpp::IntegerVector cur_indices(rIt, rIt+n);
    for (auto& i : cur_indices) { ++i; }
    slices.vec[0]=first;
    slices.vec[1]=last-first;
    
    Rcpp::Function indexed_realizer(beachenv["realizeByIndexRange"]);
    V tmp_store=indexed_realizer(original, cur_indices, slices.vec);
    std::copy(tmp_store.begin(), tmp_store.end(), out);
    return;
}

template<typename T, class V>
template<class Iter>
void unknown_reader<T, V>::get_cols(Rcpp::IntegerVector::iterator cIt, size_t n, Iter out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);

    // Need to make a copy of the indices (1-indexed) to pass to the function.
    // Don't use row_indices.vec, avoid bugs upon interaction with update_storage().
    Rcpp::IntegerVector cur_indices(cIt, cIt+n);
    for (auto& i : cur_indices) { ++i; }
    slices.vec[0]=first;
    slices.vec[1]=last-first;

    Rcpp::Function indexed_realizer(beachenv["realizeByRangeIndex"]);
    V tmp_store=indexed_realizer(original, slices.vec, cur_indices);
    std::copy(tmp_store.begin(), tmp_store.end(), out);
    return;
}

}

#endif
