#ifndef BEACHMAT_EXTERNAL_WRITER_H
#define BEACHMAT_EXTERNAL_WRITER_H

#include "Rcpp.h"

#include "../utils/utils.h"
#include "../utils/dim_checker.h"
#include "../utils/external.h"

#include <algorithm>

namespace beachmat {

/*************************
 * Base class definition *
 *************************/

template<typename T, class V>
class external_writer_base : public dim_checker {
protected:
    std::string cls, pkg;
    external_ptr ex; 

    void (*store) (void *, size_t, size_t, T*)=NULL;
    void (*load) (void *, size_t, size_t, T*)=NULL;
    SEXP (*report) (void *)=NULL;

public:
    external_writer_base(size_t nr, size_t nc, const std::string& Pkg, const std::string& Class) :
            dim_checker(nr, nc), cls(Class), pkg(Pkg), ex(nr, nc, Pkg, Class, get_type()) {

        auto type=get_type();

        // Define all remaining function pointers.
        auto store_name=get_external_name(cls, type, "output", "set");
        store=reinterpret_cast<void (*)(void *, size_t, size_t, T*)>(R_GetCCallable(pkg.c_str(), store_name.c_str()));

        auto load_name=get_external_name(cls, type, "output", "get");
        load=reinterpret_cast<void (*)(void *, size_t, size_t, T*)>(R_GetCCallable(pkg.c_str(), load_name.c_str()));

        auto report_name=get_external_name(cls, type, "output", "yield");
        report=reinterpret_cast<SEXP (*)(void *)>(R_GetCCallable(pkg.c_str(), report_name.c_str()));
        return;
    }
    ~external_writer_base() = default;
    external_writer_base(const external_writer_base&) = default;
    external_writer_base& operator=(const external_writer_base&) = default;
    external_writer_base(external_writer_base&&) = default;
    external_writer_base& operator=(external_writer_base&&) = default;

    // Setters:
    void set(size_t r, size_t c, T val) {
        store(ex.get(), r, c, &val);
        return;
    }

    // Getters:
    T get(size_t r, size_t c) {
        T tmp;
        load(ex.get(), r, c, &tmp);
        return tmp;
    }
    
    // Other:
    Rcpp::RObject yield() {
        return report(ex.get());
    }

    std::string get_class() const {
        return cls;
    }

    std::string get_package() const {
        return pkg;
    }

    static std::string get_type();
};

/*******************************
 *** Class with solo getters ***
 *******************************/

template<typename T, class V>
class external_writer : public external_writer_base<T, V> {
private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef typename V::iterator RcppValIt;

    void (*store_col) (void *, size_t, RcppValIt*, size_t, size_t);
    void (*store_row) (void *, size_t, RcppValIt*, size_t, size_t);

    void (*store_col_indexed) (void *, size_t, size_t, RcppIntIt*, RcppValIt*);
    void (*store_row_indexed) (void *, size_t, size_t, RcppIntIt*, RcppValIt*);

    void (*load_col) (void *, size_t, RcppValIt*, size_t, size_t);
    void (*load_row) (void *, size_t, RcppValIt*, size_t, size_t);

public:    
    external_writer(size_t nr, size_t nc, const std::string& Pkg, const std::string& Class) :
            external_writer_base<T, V>(nr, nc, Pkg, Class) { 

        auto Type=this->get_type();

        // Getting all required functions from the corresponding shared library.
        auto store_col_name=get_external_name(Class, Type, "output", "setCol");
        store_col=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), store_col_name.c_str()));
    
        auto store_row_name=get_external_name(Class, Type, "output", "setRow");
        store_row=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), store_row_name.c_str()));

        auto store_col_indexed_name=get_external_name(Class, Type, "output", "setColIndexed");
        store_col_indexed=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppValIt*)>(R_GetCCallable(Pkg.c_str(), store_col_indexed_name.c_str()));
    
        auto store_row_indexed_name=get_external_name(Class, Type, "output", "setRowIndexed");
        store_row_indexed=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppValIt*)>(R_GetCCallable(Pkg.c_str(), store_row_indexed_name.c_str()));

        auto load_col_name=get_external_name(Class, Type, "output", "getCol");
        load_col=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), load_col_name.c_str()));
    
        auto load_row_name=get_external_name(Class, Type, "output", "getRow");
        load_row=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), load_row_name.c_str()));

        return;
    }

    ~external_writer() = default;
    external_writer(const external_writer&) = default;
    external_writer& operator=(const external_writer&) = default;
    external_writer(external_writer&&) = default;
    external_writer& operator=(external_writer&&) = default;

    void set_row(size_t r, RcppValIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        store_row(this->ex.get(), r, &out, first, last);
        return;
    }

    void set_col(size_t c, RcppValIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        store_col(this->ex.get(), c, &out, first, last);
        return;
    }

    void set_col_indexed(size_t c, size_t n, RcppIntIt idx, RcppValIt in) {
        this->check_colargs(c);
        store_col_indexed(this->ex.get(), c, n, &idx, &in);
        return;
    }

    void set_row_indexed(size_t c, size_t n, RcppIntIt idx, RcppValIt in) {
        this->check_rowargs(c);
        store_row_indexed(this->ex.get(), c, n, &idx, &in);
        return;
    }

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
};

/******************************
 *** Class with LIN getters ***
 ******************************/

template<typename T, class V>
class external_lin_writer : public external_writer_base<T, V> {
private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef Rcpp::NumericVector::iterator RcppNumIt;

    void (*store_col_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*store_row_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*store_col_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
    void (*store_row_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);

    void (*store_col_indexed_int) (void *, size_t, size_t, RcppIntIt*, RcppIntIt*);
    void (*store_row_indexed_int) (void *, size_t, size_t, RcppIntIt*, RcppIntIt*);
    void (*store_col_indexed_dbl) (void *, size_t, size_t, RcppIntIt*, RcppNumIt*);
    void (*store_row_indexed_dbl) (void *, size_t, size_t, RcppIntIt*, RcppNumIt*);

    void (*load_col_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_row_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_col_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
    void (*load_row_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);

public:    
    external_lin_writer(size_t nr, size_t nc, const std::string& Pkg, const std::string& Class) :
            external_writer_base<T, V>(nr, nc, Pkg, Class) { 

        auto Type=this->get_type();

        // Getting all required functions from the corresponding shared library.
        auto store_col2int_name=get_external_name(Class, Type, "output", "setCol", "integer");
        store_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), store_col2int_name.c_str()));

        auto store_row2int_name=get_external_name(Class, Type, "output", "setRow", "integer");
        store_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), store_row2int_name.c_str()));

        auto store_col2dbl_name=get_external_name(Class, Type, "output", "setCol", "numeric");
        store_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), store_col2dbl_name.c_str()));

        auto store_row2dbl_name=get_external_name(Class, Type, "output", "setRow", "numeric");
        store_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), store_row2dbl_name.c_str()));

        auto store_col2int_indexed_name=get_external_name(Class, Type, "output", "setColIndexed", "integer");
        store_col_indexed_int=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppIntIt*)>(R_GetCCallable(Pkg.c_str(), store_col2int_indexed_name.c_str()));

        auto store_row2int_indexed_name=get_external_name(Class, Type, "output", "setRowIndexed", "integer");
        store_row_indexed_int=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppIntIt*)>(R_GetCCallable(Pkg.c_str(), store_row2int_indexed_name.c_str()));

        auto store_col2dbl_indexed_name=get_external_name(Class, Type, "output", "setColIndexed", "numeric");
        store_col_indexed_dbl=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppNumIt*)>(R_GetCCallable(Pkg.c_str(), store_col2dbl_indexed_name.c_str()));

        auto store_row2dbl_indexed_name=get_external_name(Class, Type, "output", "setRowIndexed", "numeric");
        store_row_indexed_dbl=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppNumIt*)>(R_GetCCallable(Pkg.c_str(), store_row2dbl_indexed_name.c_str()));

        auto load_col2int_name=get_external_name(Class, Type, "output", "getCol", "integer");
        load_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), load_col2int_name.c_str()));

        auto load_row2int_name=get_external_name(Class, Type, "output", "getRow", "integer");
        load_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), load_row2int_name.c_str()));

        auto load_col2dbl_name=get_external_name(Class, Type, "output", "getCol", "numeric");
        load_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), load_col2dbl_name.c_str()));

        auto load_row2dbl_name=get_external_name(Class, Type, "output", "getRow", "numeric");
        load_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(Pkg.c_str(), load_row2dbl_name.c_str()));

        return;
    }

    ~external_lin_writer() = default;
    external_lin_writer(const external_lin_writer&) = default;
    external_lin_writer& operator=(const external_lin_writer&) = default;
    external_lin_writer(external_lin_writer&&) = default;
    external_lin_writer& operator=(external_lin_writer&&) = default;

    // Basic setters
    void set_row(size_t r, RcppIntIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        store_row_int(this->ex.get(), r, &out, first, last);
        return;
    }

    void set_row(size_t r, RcppNumIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        store_row_dbl(this->ex.get(), r, &out, first, last);
        return;
    }

    void set_col(size_t c, RcppIntIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        store_col_int(this->ex.get(), c, &out, first, last);
        return;
    }

    void set_col(size_t c, RcppNumIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        store_col_dbl(this->ex.get(), c, &out, first, last);
        return;
    }

    // Indexed setters
    void set_col_indexed(size_t c, size_t n, RcppIntIt idx, RcppIntIt in) {
        this->check_colargs(c);
        store_col_indexed_int(this->ex.get(), c, n, &idx, &in);
        return;
    }

    void set_col_indexed(size_t c, size_t n, RcppIntIt idx, RcppNumIt in) {
        this->check_colargs(c);
        store_col_indexed_dbl(this->ex.get(), c, n, &idx, &in);
        return;
    }

    void set_row_indexed(size_t c, size_t n, RcppIntIt idx, RcppIntIt in) {
        this->check_rowargs(c);
        store_row_indexed_int(this->ex.get(), c, n, &idx, &in);
        return;
    }
    
    void set_row_indexed(size_t c, size_t n, RcppIntIt idx, RcppNumIt in) {
        this->check_rowargs(c);
        store_row_indexed_dbl(this->ex.get(), c, n, &idx, &in);
        return;
    }

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
};

}

#endif
