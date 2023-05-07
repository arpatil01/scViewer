#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

/**
 * @file lin_matrix.h
 *
 * Class definitions for the logical, integer or numeric (LIN) matrix.
 */

#include "Rcpp.h"
#include "ordinary_reader.h"
#include "Csparse_reader.h"
#include "utils.h"

#include <memory>
#include <algorithm>

namespace beachmat {

/**
 * @brief Virtual base class for a logical, integer or numeric (i.e., double-precision) matrix,
 *
 * This provides methods to extract rows or columns in dense form.
 * We suggest using the `read_lin_block()` function to construct instances of this class.
 */
class lin_matrix {
public:
    /**
     * Trivial constructor.
     */
    lin_matrix() {}

    virtual ~lin_matrix() = default;
    lin_matrix(const lin_matrix&) = default;
    lin_matrix& operator=(const lin_matrix&) = default;
    lin_matrix(lin_matrix&&) = default;
    lin_matrix& operator=(lin_matrix&&) = default;

    /**
     * Extract values from a column as an array of integers, restricted to a contiguous subset of rows.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of one-past-the-last row of interest.
     *
     * @return A pointer is returned to the values of `c` as integers, starting at the `first` element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value compares equal to `work`.
     */
    virtual const int* get_col(size_t c, int* work, size_t first, size_t last) = 0;

    /**
     * Extract values from a row as an array of integers, restricted to a contiguous subset of columns.
     *
     * @param r Index of the row of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of one-past-the-last column of interest.
     *
     * @return A pointer is returned to the values of `r` as integers, starting at the `first` element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    virtual const int* get_row(size_t r, int* work, size_t first, size_t last) = 0;

    /**
     * Extract values from a column as an array of doubles, restricted to a contiguous subset of rows.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of one-past-the-last row of interest.
     *
     * @return A pointer is returned to the values of `c` as doubles, starting at the `first` element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value compares equal to `work`.
     */
    virtual const double* get_col(size_t c, double* work, size_t first, size_t last) = 0;

    /**
     * Extract values from a row as an array of doubles, restricted to a contiguous subset of columns.
     *
     * @param r Index of the row of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of one-past-the-last column of interest.
     *
     * @return A pointer is returned to the values of `r` as doubles, starting at the `first` element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    virtual const double* get_row(size_t r, double* work, size_t first, size_t last) = 0;

    /**
     * Extract values from a column as an array of integers.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `c` as integers, starting at the first element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value compares equal to `work`.
     */
    const int* get_col(size_t c, int* work) {
        return get_col(c, work, 0, nrow);
    }

    /**
     * Extract values from a row as an array of integers.
     *
     * @param r Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `r` as integers, starting at the first element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    const int* get_row(size_t r, int* work) {
        return get_row(r, work, 0, ncol);        
    }

    /**
     * Extract values from a column as an array of doubles.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `c` as doubles, starting at the first element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value compares equal to `work`.
     */
    const double* get_col(size_t c, double* work) {
        return get_col(c, work, 0, nrow);
    }

    /**
     * Extract values from a row as an array of doubles.
     *
     * @param r Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `r` as doubles, starting at the first element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    const double* get_row(size_t r, double* work) {
        return get_row(r, work, 0, ncol);        
    }

    /**
     * Get the number of rows in the matrix.
     */
    size_t get_nrow() const { return nrow; }

    /**
     * Get the number of columns in the matrix.
     */
    size_t get_ncol() const { return ncol; }

    /**
     * Is the matrix sparse?
     */
    virtual bool is_sparse() const { return false; }

    /**
     * Clone the current object, returning a pointer to a copy.
     */
    std::unique_ptr<lin_matrix> clone() const {
        return std::unique_ptr<lin_matrix>(this->clone_internal());
    }

protected:
    size_t nrow=0, ncol=0;

    virtual lin_matrix* clone_internal() const = 0;
};

/**
 * @brief Virtual base class for a sparse logical, integer or numeric (i.e., double-precision) matrix.
 *
 * This provides methods to extract rows or columns in dense or sparse form.
 * We suggest using the `read_lin_sparse_block()` function to construct instances of this class.
 */
class lin_sparse_matrix : public lin_matrix {
public:
    /**
     * Trivial constructor.
     */
    lin_sparse_matrix() {}

    ~lin_sparse_matrix() = default;
    lin_sparse_matrix(const lin_sparse_matrix&) = default;
    lin_sparse_matrix& operator=(const lin_sparse_matrix&) = default;
    lin_sparse_matrix(lin_sparse_matrix&&) = default;
    lin_sparse_matrix& operator=(lin_sparse_matrix&&) = default;

    /**
     * Extract all non-zero elements in a row, restricted to a contiguous subset of columns.
     * Values are returned as integers.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of one-past-the-last column of interest.
     *
     * @return A `sparse_index` is returned containing pointers to the workspaces,
     * containing non-zero elements in `r` with column indices in `[first, last)`.
     * A copy of non-zero values and their indices is always performed into the workspaces.
     */
    virtual sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i, size_t first, size_t last) = 0;

    /**
     * Extract all non-zero elements in a column, restricted to a contiguous subset of rows.
     * Values are returned as integers.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of one-past-the-last row of interest.
     *
     * @return A `sparse_index` is returned containing pointers to non-zero elements in `c` with row indices in `[first, last)`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value's pointers compare equal to `work_x` and `work_i`.
     */
    virtual sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) = 0;

    /**
     * Extract all non-zero elements in a row, restricted to a contiguous subset of columns.
     * Values are returned as doubles.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of one-past-the-last column of interest.
     *
     * @return A `sparse_index` is returned containing pointers to the workspaces.
     * containing non-zero elements in `r` with column indices in `[first, last)`.
     * A copy of non-zero values and their indices is always performed into the workspaces.
     */
    virtual sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i, size_t first, size_t last) = 0;

    /**
     * Extract all non-zero elements in a column, restricted to a contiguous subset of rows.
     * Values are returned as doubles.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of one-past-the-last row of interest.
     *
     * @return A `sparse_index` is returned containing pointers to non-zero elements in `c` with row indices in `[first, last)`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value's pointers compare equal to `work_x` and `work_i`.
     */
    virtual sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last) = 0;

    /**
     * Extract all non-zero elements in a column, storing values as integers.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to all non-zero elements in `c`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value's pointers compare equal to `work_x` and `work_i`.
     */
    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i) {
        return get_col(c, work_x, work_i, 0, this->nrow);
    }

    /**
     * Extract all non-zero elements in a row, storing values as integers.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to the workspaces,
     * containing all non-zero elements in `r`.
     * A copy of non-zero values and their indices is always performed into the workspaces.
     */
    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i) {
        return get_row(r, work_x, work_i, 0, this->ncol);
    }

    /**
     * Extract all non-zero elements in a column, storing values as doubles.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to all non-zero elements in `c`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed iff the return value's pointers compare equal to `work_x` and `work_i`,
     */
    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i) {
        return get_col(c, work_x, work_i, 0, this->nrow);
    }

    /**
     * Extract all non-zero elements in a row, storing values as doubles.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to the workspace,
     * containing all non-zero elements in `r`.
     * A copy of non-zero values and their indices is always performed into the workspaces.
     */
    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i) {
        return get_row(r, work_x, work_i, 0, this->ncol);
    }

    bool is_sparse() const { return true; }

    /**
     * Get the number of non-zero elements in the matrix.
     */
    virtual size_t get_nnzero() const = 0;

    /**
     * Clone the current object, returning a pointer to a copy.
     */
    std::unique_ptr<lin_sparse_matrix> clone() const {
        return std::unique_ptr<lin_sparse_matrix>(this->clone_internal());
    }
protected:
    lin_sparse_matrix* clone_internal() const = 0;
};

/**
 * @brief Logical, integer or numeric matrices in the "ordinary" R format, i.e., column-major dense arrays.
 *
 * It is unlikely that this class will be constructed directly by users;
 * most applications will use `read_lin_block()` instead.
 *
 * @tparam V The class of the `Rcpp::Vector` holding the R-level data.
 */
template <class V>
class lin_ordinary_matrix : public lin_matrix {
public:
    /**
     * Constructor from an ordinary R-level matrix.
     *
     * @param mat An ordinary R matrix.
     */
    lin_ordinary_matrix(Rcpp::RObject mat) : reader(mat) {
        this->nrow = reader.get_nrow();
        this->ncol = reader.get_ncol();
        return;
    }

    ~lin_ordinary_matrix() = default;
    lin_ordinary_matrix(const lin_ordinary_matrix&) = default;
    lin_ordinary_matrix& operator=(const lin_ordinary_matrix&) = default;
    lin_ordinary_matrix(lin_ordinary_matrix&&) = default;
    lin_ordinary_matrix& operator=(lin_ordinary_matrix&&) = default;

    const int* get_col(size_t c, int* work, size_t first, size_t last);

    const int* get_row(size_t r, int* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last);       
        return work;
    }

    const double* get_col(size_t c, double* work, size_t first, size_t last);

    const double* get_row(size_t r, double* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last); 
        return work;
    }
private:
    ordinary_reader<V> reader;

    lin_ordinary_matrix<V>* clone_internal() const {
        return new lin_ordinary_matrix<V>(*this);
    }
};

using integer_ordinary_matrix = lin_ordinary_matrix<Rcpp::IntegerVector>;

template <>
inline const int* integer_ordinary_matrix::get_col(size_t c, int* work, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
inline const double* integer_ordinary_matrix::get_col(size_t c, double* work, size_t first, size_t last) {
    auto out = reader.get_col(c, first, last);
    std::copy(out, out + last - first, work);
    return work;
}

using logical_ordinary_matrix = lin_ordinary_matrix<Rcpp::LogicalVector>;

template <>
inline const int* logical_ordinary_matrix::get_col(size_t c, int* work, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
inline const double* logical_ordinary_matrix::get_col(size_t c, double* work, size_t first, size_t last) {
    auto out = reader.get_col(c, first, last);
    std::copy(out, out + last - first, work);
    return work;
}

using double_ordinary_matrix = lin_ordinary_matrix<Rcpp::NumericVector>;

template <>
inline const int* double_ordinary_matrix::get_col(size_t c, int* work, size_t first, size_t last) {
    auto out = reader.get_col(c, first, last);
    std::copy(out, out + last - first, work);
    return work;
}

template <>
inline const double* double_ordinary_matrix::get_col(size_t c, double* work, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

/**
 * @brief Sparse logical or numeric matrices in the `lgCMatrix` or `dgCMatrix` format, respectively, from the **Matrix** package.
 *
 * It is unlikely that this class will be constructed directly by users;
 * most applications will use `read_lin_block()` or `read_lin_sparse_block()` instead.
 *
 * @tparam V The class of the `Rcpp::Vector` holding the R-level data for non-zero values.
 */
template <class V, typename TIT>
class gCMatrix : public lin_sparse_matrix {
public:
    /**
     * Constructor from a `*gCMatrix`.
     *
     * @param mat A S4 object of the `dgCMatrix` or `lgCMatrix` class.
     */
    gCMatrix(Rcpp::RObject mat) : reader(mat) {
        this->nrow = reader.get_nrow();
        this->ncol = reader.get_ncol();
        return;
    }
   
    ~gCMatrix() = default;
    gCMatrix(const gCMatrix&) = default;
    gCMatrix& operator=(const gCMatrix&) = default;
    gCMatrix(gCMatrix&&) = default;
    gCMatrix& operator=(gCMatrix&&) = default;

    const int* get_col(size_t c, int* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;        
    }

    const int* get_row(size_t r, int* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last);       
        return work;
    }

    const double* get_col(size_t c, double* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;
    }

    const double* get_row(size_t r, double* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last); 
        return work;
    }
    
    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const int*>(r, work_x, work_i, first, last);
    }

    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const double*>(r, work_x, work_i, first, last);
    }

    size_t get_nnzero () const {
        return reader.get_nnzero();
    }
private:
    gCMatrix_reader<V, TIT> reader;

    gCMatrix<V, TIT>* clone_internal() const {
        return new gCMatrix<V, TIT>(*this);
    }
};

using lgCMatrix = gCMatrix<Rcpp::LogicalVector, const int*>;

template <>
inline sparse_index<const int*, int> lgCMatrix::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
inline sparse_index<const double*, int> lgCMatrix::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last) {
    return transplant<const double*>(reader.get_col(c, first, last), work_x, work_i);
}

using dgCMatrix = gCMatrix<Rcpp::NumericVector, const double*>;

template <>
inline sparse_index<const int*, int> dgCMatrix::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return transplant<const int*>(reader.get_col(c, first, last), work_x, work_i);
}

template <>
inline sparse_index<const double*, int> dgCMatrix::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

/**
 * @brief Sparse integer, logical or numeric matrices in the `SparseArraySeed` format from the **DelayedArray** package.
 *
 * It is unlikely that this class will be constructed directly by users;
 * most applications will use `read_lin_block()` or `read_lin_sparse_block()` instead.
 *
 * @tparam V The class of the `Rcpp::Vector` holding the R-level data for non-zero values.
 */
template <class V, typename TIT>
class lin_SparseArraySeed : public lin_sparse_matrix {
public:
    /**
     * Constructor from a `SparseArraySeed`.
     *
     * @param mat A S4 object of the `SparseArraySeed` class.
     */
    lin_SparseArraySeed(Rcpp::RObject mat) : reader(mat) {
        this->nrow = reader.get_nrow();
        this->ncol = reader.get_ncol();
        return;
    }

    ~lin_SparseArraySeed() = default;
    lin_SparseArraySeed(const lin_SparseArraySeed&) = default;
    lin_SparseArraySeed& operator=(const lin_SparseArraySeed&) = default;
    lin_SparseArraySeed(lin_SparseArraySeed&&) = default;
    lin_SparseArraySeed& operator=(lin_SparseArraySeed&&) = default;

    const int* get_col(size_t c, int* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;        
    }

    const int* get_row(size_t r, int* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last);       
        return work;
    }

    const double* get_col(size_t c, double* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;
    }

    const double* get_row(size_t r, double* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last); 
        return work;
    }

    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const int*>(r, work_x, work_i, first, last);
    }

    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const double*>(r, work_x, work_i, first, last);
    }

    size_t get_nnzero () const {
        return reader.get_nnzero();
    }
private:
    SparseArraySeed_reader<V, TIT> reader;

    lin_SparseArraySeed<V, TIT>* clone_internal() const {
        return new lin_SparseArraySeed<V, TIT>(*this);
    }
};

using integer_SparseArraySeed = lin_SparseArraySeed<Rcpp::IntegerVector, const int*>;

template <>
inline sparse_index<const int*, int> integer_SparseArraySeed::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
inline sparse_index<const double*, int> integer_SparseArraySeed::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last)
{
    return transplant<const double*>(reader.get_col(c, first, last), work_x, work_i);
}

using logical_SparseArraySeed = lin_SparseArraySeed<Rcpp::LogicalVector, const int*>;

template <>
inline sparse_index<const int*, int> logical_SparseArraySeed::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
inline sparse_index<const double*, int> logical_SparseArraySeed::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last)
{
    return transplant<const double*>(reader.get_col(c, first, last), work_x, work_i);
}

using double_SparseArraySeed = lin_SparseArraySeed<Rcpp::NumericVector, const double*>;

template <>
inline sparse_index<const int*, int> double_SparseArraySeed::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return transplant<const int*>(reader.get_col(c, first, last), work_x, work_i);
}

template <>
inline sparse_index<const double*, int> double_SparseArraySeed::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last)
{
    return reader.get_col(c, first, last);
}

}

#endif
