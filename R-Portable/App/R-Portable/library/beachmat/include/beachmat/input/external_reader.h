#ifndef BEACHMAT_EXTERNAL_READER_H
#define BEACHMAT_EXTERNAL_READER_H

#include "Rcpp.h"

#include "../utils/dim_checker.h"
#include "../utils/external.h"
#include "../utils/raw_structure.h"

#include <string>

namespace beachmat {

/******************************
 *** Basic class definition ***
 ******************************/

template<typename T, class V>
class external_reader_base : public dim_checker {
public:
    external_reader_base(const Rcpp::RObject& incoming) : original(incoming) {
        const auto& type=this->get_type();
        auto classinfo=get_class_package(original);
        cls=classinfo.first;
        pkg=classinfo.second;

        // Getting required functions from the corresponding shared library.
        auto load_name=get_external_name(cls, type, "input", "get");
        load=reinterpret_cast<void (*)(void *, size_t, size_t, T*)>(R_GetCCallable(pkg.c_str(), load_name.c_str()));

        ex=external_ptr(original, pkg, cls, type); // move assignment.

        // Getting the dimensions from the created object.
        auto get_dim_name=get_external_name(cls, type, "input", "dim");
        auto dimgetter=reinterpret_cast<void (*)(void*, size_t*, size_t*)>(R_GetCCallable(pkg.c_str(), get_dim_name.c_str()));
        dimgetter(ex.get(), &nrow, &ncol);
        return;
    }

    ~external_reader_base() = default;
    external_reader_base(const external_reader_base&) = default;
    external_reader_base& operator=(const external_reader_base&) = default;
    external_reader_base(external_reader_base&&) = default;
    external_reader_base& operator=(external_reader_base&&) = default;

    T get(size_t r, size_t c) {
        this->check_oneargs(r, c);
        T output;
        load(ex.get(), r, c, &output);
        return output;
    }

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

    // Miscellaneous.
    Rcpp::RObject yield() const { return original; }

    std::string get_class() const { return cls; }

    std::string get_package() const { return pkg; }

protected:
    Rcpp::RObject original;
    std::string cls, pkg;

    external_ptr ex;
    void (*load) (void *, size_t, size_t, T*);

    // Getting the type.
    static std::string get_type();
};

/*******************************
 *** Class with solo getters ***
 *******************************/

template<typename T, class V>
class external_reader : public external_reader_base<T, V> {
private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef typename V::iterator RcppValIt;

    void (*load_col) (void *, size_t, RcppValIt*, size_t, size_t);
    void (*load_row) (void *, size_t, RcppValIt*, size_t, size_t);

    void (*load_cols) (void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t);
    void (*load_rows) (void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t);

public:    
    external_reader(const Rcpp::RObject& incoming) : external_reader_base<T, V>(incoming) {
        const auto& type=this->get_type();
        const auto& cls=this->cls;
        const auto& pkg=this->pkg;
 
        // Getting all required functions from the corresponding shared library.
        auto load_col_name=get_external_name(cls, type, "input", "getCol");
        load_col=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_col_name.c_str()));
    
        auto load_row_name=get_external_name(cls, type, "input", "getRow");
        load_row=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_row_name.c_str()));
    
        auto load_cols_name=get_external_name(cls, type, "input", "getCols");
        load_cols=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_cols_name.c_str()));
    
        auto load_rows_name=get_external_name(cls, type, "input", "getRows");
        load_rows=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_rows_name.c_str()));

        return;
    }

    ~external_reader() = default;
    external_reader(const external_reader&) = default;
    external_reader& operator=(const external_reader&) = default;
    external_reader(external_reader&&) = default;
    external_reader& operator=(external_reader&&) = default;

    void get_row(size_t r, RcppValIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        load_row(this->ex.get(), r, &out, first, last);
        return;
    }

    void get_col(size_t c, RcppValIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        load_col(this->ex.get(), c, &out, first, last);
        return;
    }

    void get_rows(RcppIntIt rIt, size_t n, RcppValIt out, size_t first, size_t last) {
        this->check_rowargs(0, first, last);
        this->check_row_indices(rIt, n);
        load_rows(this->ex.get(), &rIt, n, &out, first, last);
        return;
    }

    void get_cols(RcppIntIt cIt, size_t n, RcppValIt out, size_t first, size_t last) {
        this->check_colargs(0, first, last);
        this->check_col_indices(cIt, n);
        load_cols(this->ex.get(), &cIt, n, &out, first, last);
        return;
    }
};

/******************************
 *** Class with LIN getters ***
 ******************************/

template<typename T, class V>
class external_lin_reader : public external_reader_base<T, V> {
private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef Rcpp::NumericVector::iterator RcppNumIt;

    void (*load_col_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_row_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_col_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
    void (*load_row_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);

    void (*load_cols_int) (void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t);
    void (*load_rows_int) (void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t);
    void (*load_cols_dbl) (void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t);
    void (*load_rows_dbl) (void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t);
public:    
    external_lin_reader(const Rcpp::RObject& incoming) : external_reader_base<T, V>(incoming) {
        const auto& type=this->get_type();
        const auto& cls=this->cls;
        const auto& pkg=this->pkg;

        // Getting all required functions from the corresponding shared library.
        auto load_col2int_name=get_external_name(cls, type, "input", "getCol", "integer");
        load_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_col2int_name.c_str()));

        auto load_row2int_name=get_external_name(cls, type, "input", "getRow", "integer");
        load_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_row2int_name.c_str()));

        auto load_col2dbl_name=get_external_name(cls, type, "input", "getCol", "numeric");
        load_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_col2dbl_name.c_str()));

        auto load_row2dbl_name=get_external_name(cls, type, "input", "getRow", "numeric");
        load_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_row2dbl_name.c_str()));

        auto load_cols2int_name=get_external_name(cls, type, "input", "getCols", "integer");
        load_cols_int=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_cols2int_name.c_str()));

        auto load_rows2int_name=get_external_name(cls, type, "input", "getRows", "integer");
        load_rows_int=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_rows2int_name.c_str()));

        auto load_cols2dbl_name=get_external_name(cls, type, "input", "getCols", "numeric");
        load_cols_dbl=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_cols2dbl_name.c_str()));

        auto load_rows2dbl_name=get_external_name(cls, type, "input", "getRows", "numeric");
        load_rows_dbl=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg.c_str(), load_rows2dbl_name.c_str()));

        return;
    }

    ~external_lin_reader() = default;
    external_lin_reader(const external_lin_reader&) = default;
    external_lin_reader& operator=(const external_lin_reader&) = default;
    external_lin_reader(external_lin_reader&&) = default;
    external_lin_reader& operator=(external_lin_reader&&) = default;

    // Basic getters
    void get_row(size_t r, RcppIntIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        load_row_int(this->ex.get(), r, &out, first, last);
        return;
    }

    void get_row(size_t r, RcppNumIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        load_row_dbl(this->ex.get(), r, &out, first, last);
        return;
    }

    void get_col(size_t c, RcppIntIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        load_col_int(this->ex.get(), c, &out, first, last);
        return;
    }

    void get_col(size_t c, RcppNumIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        load_col_dbl(this->ex.get(), c, &out, first, last);
        return;
    }

    // Multi getters
    void get_rows(RcppIntIt rIt, size_t n, RcppIntIt out, size_t first, size_t last) {
        this->check_rowargs(0, first, last);
        this->check_row_indices(rIt, n);
        load_rows_int(this->ex.get(), &rIt, n, &out, first, last);
        return;
    }
    
    void get_rows(RcppIntIt rIt, size_t n, RcppNumIt out, size_t first, size_t last) {
        this->check_rowargs(0, first, last);
        this->check_row_indices(rIt, n);
        load_rows_dbl(this->ex.get(), &rIt, n, &out, first, last);
        return;
    }
    
    void get_cols(RcppIntIt cIt, size_t n, RcppIntIt out, size_t first, size_t last) {
        this->check_colargs(0, first, last);
        this->check_col_indices(cIt, n);
        load_cols_int(this->ex.get(), &cIt, n, &out, first, last);
        return;
    }
    
    void get_cols(RcppIntIt cIt, size_t n, RcppNumIt out, size_t first, size_t last) {
        this->check_colargs(0, first, last);
        this->check_col_indices(cIt, n);
        load_cols_dbl(this->ex.get(), &cIt, n, &out, first, last);
        return;
    }
};

}

#endif
