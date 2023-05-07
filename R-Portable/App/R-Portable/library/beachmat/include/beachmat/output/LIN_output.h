#ifndef BEACHMAT_LIN_OUTPUT_H
#define BEACHMAT_LIN_OUTPUT_H

#include "Rcpp.h"

#include "simple_writer.h"
#include "Csparse_writer.h"
#include "external_writer.h"
#include "output_param.h"
#include "../utils/utils.h"

#include <memory>

namespace beachmat { 

/************************************************************************
 * Virtual base class for LIN (logical/integer/numeric) output matrices.
 ************************************************************************/

template<typename T, class V>
class lin_output {
public:
    lin_output() = default;
    virtual ~lin_output() = default;
    lin_output(const lin_output&) = default;
    lin_output& operator=(const lin_output&) = default;
    lin_output(lin_output&&) = default;
    lin_output& operator=(lin_output&&) = default;
    
    // Getters:
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;

    void get_row(size_t, Rcpp::IntegerVector::iterator);
    void get_row(size_t, Rcpp::NumericVector::iterator);

    virtual void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_col(size_t, Rcpp::IntegerVector::iterator);
    void get_col(size_t, Rcpp::NumericVector::iterator);

    virtual void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    virtual T get(size_t, size_t)=0;

    // Setters:
    void set_row(size_t, Rcpp::IntegerVector::iterator);
    void set_row(size_t, Rcpp::NumericVector::iterator);

    virtual void set_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void set_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void set_col(size_t, Rcpp::IntegerVector::iterator);
    void set_col(size_t, Rcpp::NumericVector::iterator);

    virtual void set_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void set_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    virtual void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::IntegerVector::iterator)=0;
    virtual void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::NumericVector::iterator)=0;

    virtual void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::IntegerVector::iterator)=0;
    virtual void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::NumericVector::iterator)=0;

    virtual void set(size_t, size_t, T)=0;

    // Other methods:
    virtual Rcpp::RObject yield()=0;

    virtual std::unique_ptr<lin_output<T, V> > clone() const=0;

    virtual std::string get_class() const=0;

    virtual std::string get_package() const=0;

    // Useful typedefs
    typedef V vector;

    typedef T type;
private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed.
};

/* General output */

template<typename T, class V, class WTR>
class general_lin_output : public lin_output<T, V> {
public:
    general_lin_output(size_t, size_t);
    ~general_lin_output() = default;
    general_lin_output(const general_lin_output&) = default;
    general_lin_output& operator=(const general_lin_output&) = default;
    general_lin_output(general_lin_output&&) = default;
    general_lin_output& operator=(general_lin_output&&) = default;

    // Getters:
    size_t get_nrow() const;
    size_t get_ncol() const;

    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    // Setters:
    void set_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void set_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void set_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void set_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::IntegerVector::iterator);
    void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::NumericVector::iterator);

    void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::IntegerVector::iterator);
    void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::NumericVector::iterator);

    void set(size_t, size_t, T);

    // Other:
    Rcpp::RObject yield();

    std::unique_ptr<lin_output<T, V> > clone() const;

    std::string get_class() const { return writer.get_class(); }

    std::string get_package() const { return writer.get_package(); }
protected:
    general_lin_output(WTR&& w) : writer(w) {}
    WTR writer;
};

/* Simple LIN output */

template<typename T, class V>
using simple_lin_output=general_lin_output<T, V, simple_writer<T, V> >;

/* Sparse LIN output */

template<typename T, class V>
using sparse_lin_output=general_lin_output<T, V, Csparse_writer<T, V> >;

/* External LIN output */

template<typename T, class V>
class external_lin_output : public general_lin_output<T, V, external_lin_writer<T, V> > {
public:
    external_lin_output(size_t nr, size_t nc, const std::string& pkg, const std::string& cls) :
        general_lin_output<T, V, external_lin_writer<T, V> >(external_lin_writer<T, V>(nr, nc, pkg, cls)) {}
    ~external_lin_output() = default;
    external_lin_output(const external_lin_output&) = default;
    external_lin_output& operator=(const external_lin_output&) = default;
    external_lin_output(external_lin_output&&) = default;
    external_lin_output& operator=(external_lin_output&&) = default;
};

}

#include "LIN_methods.h"

#endif

