#ifndef BEACHMAT_DIM_CHECKER_H
#define BEACHMAT_DIM_CHECKER_H

#include "Rcpp.h"
#include <stdexcept>

/**
 * @file dim_checker.h
 *
 * Internal class definition for validating requested indices against the data dimensionality.
 */

namespace beachmat {

/**
 * @brief Base virtual class implementing dimensionality checks.
 *
 * @note This is an internal class and should not be constructed directly by **beachmat** users.
 */
class dim_checker {
public:
    dim_checker() {}

    virtual ~dim_checker() = default;
    dim_checker(const dim_checker&) = default;
    dim_checker& operator=(const dim_checker&) = default;
    dim_checker(dim_checker&&) = default;
    dim_checker& operator=(dim_checker&&) = default;

    /**
     * Check index against the dimension length and report an error if the former is out of range.
     *
     * @param i Requested index along a dimension of interest.
     * @param dim Total length of the dimension.
     * @param msg Name of the dimension, e.g., `"row"`.
     *
     * @return An error is raised if `i >= dim`.
     */
    static void check_dimension(size_t i, size_t dim, const std::string& msg) {
        if (i >= dim) {
            throw std::runtime_error(msg + " index out of range");
        }
        return;
    }

    /** 
     * Check requested subset indices against the dimension and report errors if the former are invalid.
     *
     * @param first Requested start index along a dimension of interest.
     * @param last Requested one-past-the-end index along the dimension of interest.
     * @param dim Total length of the dimension.
     * @param msg Name of the dimension, e.g., `"row"`.
     *
     * @return An error is raised if `last < first` or if `last > dim`.
     */
    static void check_subset(size_t first, size_t last, size_t dim, const std::string& msg) {
         if (last < first) {
            throw std::runtime_error(msg + " start index is greater than " + msg + " end index");
         } else if (last > dim) {
             throw std::runtime_error(msg + " end index out of range");
         }
         
         return;    
    }

    /** 
     * Get the current number of rows.
     */
    size_t get_nrow() const { return nrow; }

    /** 
     * Get the current number of columns.
     */
    size_t get_ncol() const { return ncol; }
protected:
    size_t nrow=0, ncol=0;

    /**
     * Fill the dimensions of this object.
     *
     * @param dims An `Rcpp::IntegerVector` of length 2 or something that can be coerced to one.
     * 
     * @return `nrow` and `ncol` are set to the first and second element of `dims`, respectively.
     */
    void fill_dims(Rcpp::RObject dims) {
        if (dims.sexp_type()!=INTSXP) {
            throw std::runtime_error("matrix dimensions should be an integer vector");
        }

        Rcpp::IntegerVector d(dims);
        if (d.size()!=2) {
            throw std::runtime_error("matrix dimensions should be of length 2");
        }

        if (d[0]<0 || d[1]<0) {
            throw std::runtime_error("dimensions should be non-negative");
        }
        nrow=d[0];
        ncol=d[1];
        return;
    }

    /**
     * Check that the requested row index is compatible with the stored dimensions.
     *
     * @param r A requested row index.
     *
     * @return An error is raised for invalid `r`, see `check_dimension()`.
     */
    void check_rowargs(size_t r) const {
        dim_checker::check_dimension(r, nrow, "row");
        return;
    }

    /**
     * Check that the requested row index and column subsets are compatible with the stored dimensions.
     *
     * @param r A requested row index.
     * @param first Index of the first column of interest.
     * @param last Index of one-past-the-last column of interest.
     *
     * @return An error is raised for invalid indices, see `check_dimension()` and `check_subset()`.
     */
    void check_rowargs(size_t r, size_t first, size_t last) const {
        check_rowargs(r);
        dim_checker::check_subset(first, last, ncol, "column");
        return;
    }

    /**
     * Check that the requested column index is compatible with the stored dimensions.
     *
     * @param c A requested column index.
     *
     * @return An error is raised for invalid `c`, see `check_dimension()`.
     */
    void check_colargs(size_t c) const {
        dim_checker::check_dimension(c, ncol, "column");
        return;
    }

    /**
     * Check that the requested column index and row subsets are compatible with the stored dimensions.
     *
     * @param c A requested column index.
     * @param first Index of the first row of interest.
     * @param last Index of one-past-the-last row of interest.
     *
     * @return An error is raised for invalid indices, see `check_dimension()` and `check_subset()`.
     */
    void check_colargs(size_t c, size_t first, size_t last) const {
        check_colargs(c);
        dim_checker::check_subset(first, last, nrow, "row");
        return;
    }
};

}

#endif
