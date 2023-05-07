#ifndef BEACHMAT_CHARACTER_MATRIX_H
#define BEACHMAT_CHARACTER_MATRIX_H

#include "Rcpp.h"

#include "input/simple_reader.h"
#include "input/dense_reader.h"
#include "input/delayed_reader.h"
#include "input/unknown_reader.h"
#include "input/external_reader.h"

#include "output/simple_writer.h"
#include "output/output_param.h"

#include "utils/utils.h"

#include <memory>
#include <vector>

namespace beachmat { 

/*********
 * INPUT *
 *********/

/* Virtual base class for character matrices. */

class character_matrix {
public:    
    character_matrix() = default;
    virtual ~character_matrix() = default;
    character_matrix(const character_matrix&) = default;
    character_matrix& operator=(const character_matrix&) = default;
    character_matrix(character_matrix&&) = default;
    character_matrix& operator=(character_matrix&&) = default;

    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;
    
    // Basic getters.
    void get_col(size_t c, Rcpp::StringVector::iterator out) { 
        get_col(c, out, 0, get_nrow());
        return;
    }
    virtual void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_row(size_t r, Rcpp::StringVector::iterator out) { 
        get_row(r, out, 0, get_ncol());
        return;
    }
    virtual void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual Rcpp::String get(size_t, size_t)=0;

    // Specialized getters.
    virtual raw_structure<Rcpp::StringVector> set_up_raw () const=0;

    void get_col_raw(size_t c, raw_structure<Rcpp::StringVector>& in) {
        get_col_raw(c, in, 0, get_nrow());
        return;
    }
    virtual void get_col_raw(size_t, raw_structure<Rcpp::StringVector>&, size_t, size_t)=0;

    void get_row_raw(size_t r, raw_structure<Rcpp::StringVector>& in) {
        get_row_raw(r, in, 0, get_ncol());
        return;
    }
    virtual void get_row_raw(size_t, raw_structure<Rcpp::StringVector>&, size_t, size_t)=0;

    virtual std::string col_raw_type() const = 0;

    virtual std::string row_raw_type() const = 0;

    // Multi getters.
    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out) {
        get_cols(it, n, out, 0, get_nrow());
        return;
    }
    virtual void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out) {
        get_rows(it, n, out, 0, get_ncol());
        return;
    }
    virtual void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    // Other methods.
    virtual std::unique_ptr<character_matrix> clone() const=0;

    virtual Rcpp::RObject yield () const=0;

    virtual std::string get_class() const=0;

    virtual std::string get_package() const=0;

    // Useful typedefs
    typedef Rcpp::StringVector vector;

    typedef Rcpp::String type;
};

std::unique_ptr<character_matrix> create_character_matrix_internal(const Rcpp::RObject&, bool); 

/* Advanced character matrix template */

template<class RDR>
class general_character_matrix : public character_matrix {
public:    
    general_character_matrix(const Rcpp::RObject& incoming) : reader(incoming) {}
    ~general_character_matrix() = default;
    general_character_matrix(const general_character_matrix&) = default;
    general_character_matrix& operator=(const general_character_matrix&) = default;
    general_character_matrix(general_character_matrix&&) = default;
    general_character_matrix& operator=(general_character_matrix&&) = default;
  
    size_t get_nrow() const { return reader.get_nrow(); }
    size_t get_ncol() const { return reader.get_ncol(); }

    // Basic getters.
    using character_matrix::get_row;
    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        reader.get_row(r, out, first, last);
        return;
    }

    using character_matrix::get_col;
    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        reader.get_col(c, out, first, last); 
        return;
    }

    Rcpp::String get(size_t r, size_t c) { return reader.get(r, c); }

    // Specialized getters.
    raw_structure<Rcpp::StringVector> set_up_raw() const {
        return reader.set_up_raw();
    }

    using character_matrix::get_col_raw;
    void get_col_raw(size_t c, raw_structure<Rcpp::StringVector>& in, size_t first, size_t last) {
        reader.get_col_raw(c, in, first, last);
        return;
    }
    
    using character_matrix::get_row_raw;
    void get_row_raw(size_t r, raw_structure<Rcpp::StringVector>& in, size_t first, size_t last) {
        reader.get_row_raw(r, in, first, last);
        return;
    }

    virtual std::string col_raw_type() const {
        return reader.col_raw_type();
    }

    virtual std::string row_raw_type() const {
        return reader.row_raw_type();
    }

    // Multi getters.
    using character_matrix::get_rows;
    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        reader.get_rows(it, n, out, first, last);
        return;
    }

    using character_matrix::get_cols;
    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        reader.get_cols(it, n, out, first, last);
        return;
    }

    // Other methods.
    std::unique_ptr<character_matrix> clone() const { return std::unique_ptr<character_matrix>(new general_character_matrix(*this)); }

    Rcpp::RObject yield () const { return reader.yield(); }

    std::string get_class() const { return reader.get_class(); }

    std::string get_package() const { return reader.get_package(); }
protected:
    RDR reader;
};

/* Simple character matrix */

using simple_character_matrix=general_character_matrix<simple_reader<Rcpp::String, Rcpp::StringVector> >;

/* DelayedMatrix */

typedef delayed_reader<Rcpp::String, Rcpp::StringVector, character_matrix> delayed_character_reader;

template<>
inline std::unique_ptr<character_matrix> delayed_character_reader::generate_seed(Rcpp::RObject incoming) {
    return create_character_matrix_internal(incoming, false);
}

using delayed_character_matrix=general_character_matrix<delayed_character_reader>;

/* Unknown matrix type */

using unknown_character_matrix=general_character_matrix<unknown_reader<Rcpp::String, Rcpp::StringVector> >;

/* External matrix type */

template<>
inline std::string external_reader_base<Rcpp::String, Rcpp::StringVector>::get_type() { return "character"; }

using external_character_matrix=general_character_matrix<external_reader<Rcpp::String, Rcpp::StringVector> >;

/* Dispatcher */

inline std::unique_ptr<character_matrix> create_character_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) { 
        std::string ctype=get_class_name(incoming);
        if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<character_matrix>(new delayed_character_matrix(incoming));
        } else if (has_external_support("character", incoming)) {
            return std::unique_ptr<character_matrix>(new external_character_matrix(incoming));
        }
        return std::unique_ptr<character_matrix>(new unknown_character_matrix(incoming));
    } 
    quit_on_df(incoming);
    return std::unique_ptr<character_matrix>(new simple_character_matrix(incoming));
}

inline std::unique_ptr<character_matrix> create_character_matrix(const Rcpp::RObject& incoming) { 
    return create_character_matrix_internal(incoming, true);
}

template<>
inline std::unique_ptr<character_matrix> create_matrix<character_matrix>(const Rcpp::RObject& incoming) {
    return create_character_matrix(incoming);
}

/**********
 * OUTPUT *
 **********/

/* Virtual base class for character matrices. */

class character_output {
public:    
    character_output() = default;
    virtual ~character_output() = default;
    character_output(const character_output&) = default;
    character_output& operator=(const character_output&) = default;
    character_output(character_output&&) = default;
    character_output& operator=(character_output&&) = default;
    
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;

    // Getters    
    void get_col(size_t c, Rcpp::StringVector::iterator out) { 
        get_col(c, out, 0, get_nrow());
        return;
    }
    virtual void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_row(size_t r, Rcpp::StringVector::iterator out) { 
        get_row(r, out, 0, get_ncol());
        return;
    }
    virtual void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual Rcpp::String get(size_t, size_t)=0;

    // Setters
    void set_col(size_t c, Rcpp::StringVector::iterator out) { 
        set_col(c, out, 0, get_nrow());
        return;
    }
    virtual void set_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void set_row(size_t r, Rcpp::StringVector::iterator out) { 
        set_row(r, out, 0, get_ncol());
        return;
    }
    virtual void set_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual void set(size_t, size_t, Rcpp::String)=0;

    virtual void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::StringVector::iterator)=0;

    virtual void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Rcpp::StringVector::iterator)=0;

    // Other stuff.
    virtual Rcpp::RObject yield()=0;

    virtual std::unique_ptr<character_output> clone() const=0;

    virtual std::string get_class() const=0;

    virtual std::string get_package() const=0;

    // Useful typedefs
    typedef Rcpp::StringVector vector;

    typedef Rcpp::String type;
};

/* General character matrix */

template<class WTR>
class general_character_output : public character_output {
public:
    general_character_output(size_t nr, size_t nc) : writer(nr, nc) {}
    ~general_character_output() = default;
    general_character_output(const general_character_output&) = default;
    general_character_output& operator=(const general_character_output&) = default;
    general_character_output(general_character_output&&) = default;
    general_character_output& operator=(general_character_output&&) = default;

    size_t get_nrow() const {
        return writer.get_nrow();
    }
    
    size_t get_ncol() const {
        return writer.get_ncol();
    }
    
    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        writer.get_row(r, out, first, last);
        return;
    }
    
    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
        writer.get_col(c, out, first, last);
        return;
    }
    
    Rcpp::String get(size_t r, size_t c) {
        return writer.get(r, c);
    }
    
    void set_row(size_t r, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
        writer.set_row(r, in, first, last);
        return;
    }
    
    void set_col(size_t c, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
        writer.set_col(c, in, first, last);
        return;
    }
    
    void set(size_t r, size_t c, Rcpp::String in) {
        writer.set(r, c, in);
        return;
    }
    
    void set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
        writer.set_col_indexed(c, N, idx, val);
        return;
    }
    
    void set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
        writer.set_row_indexed(r, N, idx, val);
        return;
    }
    
    Rcpp::RObject yield() {
        return writer.yield();
    }
    
    std::unique_ptr<character_output> clone() const {
        return std::unique_ptr<character_output>(new general_character_output<WTR>(*this));
    }

    std::string get_class() const {
        return writer.get_class();
    }

    std::string get_package() const {
        return writer.get_package();
    }
protected:
    general_character_output(WTR&& w) : writer(w) {}
    WTR writer;
};

/* Simple character matrix */

typedef general_character_output<simple_writer<Rcpp::String, Rcpp::StringVector> > simple_character_output;

/* External character matrix */

template<>
inline std::string external_writer_base<Rcpp::String, Rcpp::StringVector>::get_type() { return "character"; }

class external_character_output : public general_character_output<external_writer<Rcpp::String, Rcpp::StringVector> > {
public:
    external_character_output(size_t nr, size_t nc, const std::string& pkg, const std::string& cls) :
        general_character_output<external_writer<Rcpp::String, Rcpp::StringVector> >(
            external_writer<Rcpp::String, Rcpp::StringVector>(nr, nc, pkg, cls)) {}
    ~external_character_output() = default;
    external_character_output(const external_character_output&) = default;
    external_character_output& operator=(const external_character_output&) = default;
    external_character_output(external_character_output&&) = default;
    external_character_output& operator=(external_character_output&&) = default;

    std::unique_ptr<character_output> clone() const {
        return std::unique_ptr<character_output>(new external_character_output(*this));
    }

};

/* Dispatcher */

inline std::unique_ptr<character_output> create_character_output(int nrow, int ncol, const output_param& param) {
    if (param.is_external_available("character")) { 
        return std::unique_ptr<character_output>(new external_character_output(nrow, ncol, param.get_package(), param.get_class()));
    }

    return std::unique_ptr<character_output>(new simple_character_output(nrow, ncol));
}

template<>
inline std::unique_ptr<character_output> create_output<character_output>(int nrow, int ncol, const output_param& param) {
    return create_character_output(nrow, ncol, param);
}

}

#endif
