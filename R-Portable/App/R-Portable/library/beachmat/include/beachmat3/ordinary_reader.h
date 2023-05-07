#ifndef BEACHMAT_ORDINARY_READER_H
#define BEACHMAT_ORDINARY_READER_H

/**
 * @file ordinary_reader.h
 *
 * Internal class definition to read from ordinary matrices.
 */

#include "Rcpp.h"
#include "dim_checker.h"
#include "utils.h"

namespace beachmat {

/**
 * @brief Type-agnostic reader for row/column data from an ordinary R matrix.
 *
 * Unlike the `lin_ordinary_matrix` class template, this is type-agnostic.
 *
 * @note This is an internal class and should not be constructed directly by **beachmat** users.
 *
 * @tparam V An `Rcpp::Vector` class, used to contain the contents of the matrix.
 */
template <class V>
class ordinary_reader : public dim_checker {
public:
    ~ordinary_reader() = default;
    ordinary_reader(const ordinary_reader&) = default;
    ordinary_reader& operator=(const ordinary_reader&) = default;
    ordinary_reader(ordinary_reader&&) = default;
    ordinary_reader& operator=(ordinary_reader&&) = default;

    /** 
     * Constructor from an ordinary R matrix.
     *
     * @param An ordinary R matrix.
     */
    ordinary_reader(Rcpp::RObject input) : mat(input) {
        this->fill_dims(input.attr("dim")); // consistency between 'dim' and mat.size() is guaranteed by R.
        return;
    }

    /**
     * Return an iterator to a sequence of values from a column of the matrix,
     * possibly restricted to a subset of rows.
     *
     * @param c The index of the column to extract.
     * @param first Index of the first row of interest.
     * @param last Index of one-past-the-last row of interest.
     *
     * @return An iterator pointing to the `first` element in column `c`.
     */
    typename V::const_iterator get_col(size_t c, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        return mat.begin() + (c * this->nrow) + first;
    }

    /**
     * Extract values from a row of the matrix, possibly restricted to a subset of columns.
     *
     * @param r The index of the row to extract.
     * @param workspace A pointer to an array of values in which to store the row values.
     * This should have at least `last - first` accessible elements.
     * @param first Index of the first column of interest.
     * @param last Index of one-past-the-last column of interest.
     *
     * @tparam Iter A pointer or iterator to a writeable sequence of elements.
     *
     * @return Row values from row `r` and from columns `[first, last)` are copied to the `workspace`.
     */
    template <class Iter>
    void get_row(size_t r, Iter work, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        auto src=mat.begin() + first * (this->nrow) + r;
        for (size_t col=first; col<last; ++col, src+=(this->nrow), ++work) { (*work)=(*src); }
        return;
    }
private:
    V mat;
};

}

#endif
