#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

#include "Rcpp.h"

#include "simple_reader.h"
#include "dense_reader.h"
#include "Csparse_reader.h"
#include "delayed_reader.h"
#include "unknown_reader.h"
#include "external_reader.h"

#include "../utils/utils.h"
#include "../utils/raw_structure.h"

#include <memory>

namespace beachmat { 

/***************************************************************** 
 * Virtual base class for LIN (logical/integer/numeric) matrices. 
 *****************************************************************/

template<typename T, class V>
class lin_matrix {
public:
    lin_matrix() = default;
    virtual ~lin_matrix() = default;
    lin_matrix(const lin_matrix&) = default;
    lin_matrix& operator=(const lin_matrix&) = default;
    lin_matrix(lin_matrix&&) = default;
    lin_matrix& operator=(lin_matrix&&) = default;
    
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;

    /* We can't add a LogicalVector::iterator method because IntegerVector::iterator==LogicalVector::iterator
     * under the hood in Rcpp. The compiler then complains that overloading is not possible. Thus, for all 
     * references here to LogicalVector, we will consider the use of IntegerVector in its place.
     */

    void get_row(size_t, Rcpp::IntegerVector::iterator);
    void get_row(size_t, Rcpp::NumericVector::iterator);

    virtual void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_col(size_t, Rcpp::IntegerVector::iterator);
    void get_col(size_t, Rcpp::NumericVector::iterator);

    virtual void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    virtual T get(size_t, size_t)=0;
   
    // Specialized getters.
    virtual raw_structure<V> set_up_raw() const=0;
  
    void get_col_raw(size_t c, raw_structure<V>& in) {
        get_col_raw(c, in, 0, get_nrow());
        return;
    }

    virtual void get_col_raw(size_t, raw_structure<V>&, size_t, size_t)=0;

    void get_row_raw(size_t r, raw_structure<V>& in) {
        get_row_raw(r, in, 0, get_ncol());
        return;
    }

    virtual void get_row_raw(size_t, raw_structure<V>&, size_t, size_t)=0;

    virtual std::string col_raw_type() const = 0;

    virtual std::string row_raw_type() const = 0;

    // Multi-row/column getters.
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator);
    void get_rows(Rcpp::IntegerVector::iterator ,size_t, Rcpp::NumericVector::iterator);

    virtual void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator);

    virtual void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    // Other methods.
    virtual std::unique_ptr<lin_matrix<T, V> > clone() const=0;

    virtual Rcpp::RObject yield() const=0;

    virtual std::string get_class() const=0;

    virtual std::string get_package() const=0;

    // Useful typedefs.
    typedef V vector;

    typedef T type;
};

/* A general flavour for a LIN matrix */

template <typename T, class V, class RDR>
class general_lin_matrix : public lin_matrix<T, V> {
public:    
    general_lin_matrix(const Rcpp::RObject&);
    ~general_lin_matrix() = default;
    general_lin_matrix(const general_lin_matrix&) = default;
    general_lin_matrix& operator=(const general_lin_matrix&) = default;
    general_lin_matrix(general_lin_matrix&&) = default;
    general_lin_matrix& operator=(general_lin_matrix&&) = default;
    
    size_t get_nrow() const;
    size_t get_ncol() const;

    // Basic getters.
    using lin_matrix<T, V>::get_col;
    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    using lin_matrix<T, V>::get_row;
    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    // Specialized getters (these do nothing by default, see specializations below).
    raw_structure<V> set_up_raw() const {
        return reader.set_up_raw();
    }

    using lin_matrix<T, V>::get_col_raw;
    void get_col_raw(size_t c, raw_structure<V>& in, size_t first, size_t last) {
        reader.get_col_raw(c, in, first, last);
        return;
    }
    
    using lin_matrix<T, V>::get_row_raw;
    void get_row_raw(size_t r, raw_structure<V>& in, size_t first, size_t last) {
        reader.get_row_raw(r, in, first, last);
        return;
    }
 
    virtual std::string col_raw_type() const {
        return reader.col_raw_type();
    }

    virtual std::string row_raw_type() const {
        return reader.row_raw_type();
    }

    // Multigetters.
    using lin_matrix<T, V>::get_rows;
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    using lin_matrix<T, V>::get_cols;
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;

    Rcpp::RObject yield() const;

    std::string get_class() const { return reader.get_class(); }

    std::string get_package() const { return reader.get_package(); }
protected:
    RDR reader;
};

/* Realization of the general flavour for a simple matrix */

template <typename T, class V>
using simple_lin_matrix=general_lin_matrix<T, V, simple_reader<T, V> >;

/* Realization of the general flavour for a dense matrix */

template <typename T, class V>
using dense_lin_matrix=general_lin_matrix<T, V, dense_reader<T, V> >;

/* Realization of the general flavour for a C-sparse matrix */

template <typename T, class V>
using Csparse_lin_matrix=general_lin_matrix<T, V, Csparse_reader<T, V> >;

/* DelayedMatrix of LINs */

template <typename T, class V>
using delayed_lin_reader=delayed_reader<T, V, lin_matrix<T, V> >;

template <typename T, class V>
using delayed_lin_matrix=general_lin_matrix<T, V, delayed_lin_reader<T, V> >;

/* Unknown matrix of LINs */

template <typename T, class V>
using unknown_lin_matrix=general_lin_matrix<T, V, unknown_reader<T, V> >;

/* External matrix of LINs */

template <typename T, class V>
using external_lin_matrix=general_lin_matrix<T, V, external_lin_reader<T, V> >;

}

#include "LIN_methods.h"

#endif

