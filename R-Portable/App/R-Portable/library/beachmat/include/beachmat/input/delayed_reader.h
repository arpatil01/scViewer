#ifndef BEACHMAT_DELAYED_READER_H
#define BEACHMAT_DELAYED_READER_H

#include "Rcpp.h"

#include "unknown_reader.h"
#include "../utils/utils.h"
#include "../utils/dim_checker.h"
#include "../utils/copyable_vector.h"
#include "../utils/raw_structure.h"

#include <memory>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace beachmat {

/* The 'delayed_coord_transformer' class is intended to allow direct data access
 * from the underlying seed matrix in a DelayedMatrix class, while transforming
 * the extracted values to account for row/column subsetting or transposition.
 * This avoids the need to realize any subset of the matrix, as would be necessary
 * for more general delayed operations (see the 'delayed_reader' class below).
 */

template<typename T, class V>
class delayed_coord_transformer {
public:    
    delayed_coord_transformer() = default;

    template<class M>
    delayed_coord_transformer(M);

    template<class M>
    delayed_coord_transformer(const Rcpp::List&, const Rcpp::LogicalVector&, M);

    template<class M, class Iter>
    void get_row(M, size_t, Iter, size_t, size_t);
    
    template<class M, class Iter>
    void get_col(M, size_t, Iter, size_t, size_t);

    template<class M>
    T get(M, size_t, size_t);

    size_t get_nrow() const;
    size_t get_ncol() const;
private:
    std::vector<size_t> row_index, col_index;
    bool transposed=false, byrow=false, bycol=false;
    size_t delayed_nrow=0, delayed_ncol=0;

    static void obtain_indices(const Rcpp::RObject&, size_t, bool&, size_t&, std::vector<size_t>&);
    copyable_holder<V> tmp;

    // Various helper functions to implement the effect of the delayed subsetting.
    template<class M, class Iter>
    void reallocate_row(M, size_t, size_t, size_t, Iter out);
    template<class M, class Iter>
    void reallocate_col(M, size_t, size_t, size_t, Iter out);

    size_t old_col_first=0, old_col_last=0, min_col_index=0, max_col_index=0;
    size_t old_row_first=0, old_row_last=0, min_row_index=0, max_row_index=0;
    static void prepare_reallocation(size_t, size_t, size_t&, size_t&, size_t&, size_t&, const std::vector<size_t>&, const char*);
};

/* The 'delayed_reader' class, which wraps the coord_transformer class. */

template<typename T, class V, class base_mat>
class delayed_reader : public dim_checker { 
public:
    delayed_reader(const Rcpp::RObject& incoming) : original(incoming), seed_ptr(nullptr) {
        auto classinfo=get_class_package(incoming);
        if (classinfo.first!=get_class() || classinfo.second!=get_package()) {
            throw std::runtime_error("input matrix should be a DelayedMatrix");
        }

        // Parsing the delayed operation structure.
        const Rcpp::Environment beachenv=Rcpp::Environment::namespace_env("beachmat");
        Rcpp::Function parser(beachenv["setupDelayedMatrix"]);
        Rcpp::List parse_out=parser(incoming);
        if (parse_out.size()!=3) {
            throw std::runtime_error("output of beachmat:::setupDelayedMatrix should be a list of length 3");
        }

        /* Checking the matrix does not have value-operating delayed operations, 
         * and thus we can use native methods for extraction. Otherwise,
         * parsed_mat would be a DelayedMatrix (note that generate_seed
         * will simply default to an unknown matrix in such cases, rather 
         * than infinitely recursing through the delayed_matrix constructors).
         */
        Rcpp::RObject parsed_mat=parse_out[2];
        seed_ptr=generate_seed(parsed_mat);

        bool direct_extract=true;
        if (parsed_mat.isS4()) {
            auto parsedcls=get_class_package(parsed_mat);
            direct_extract=(parsedcls.first!=get_class() || parsedcls.second!=get_package());
        }

        if (direct_extract) { 
           transformer=delayed_coord_transformer<T, V>(parse_out[0], parse_out[1], seed_ptr.get());
        } else {
            // otherwise, treat it as an unknown matrix.
            transformer=delayed_coord_transformer<T, V>(seed_ptr.get());
        }

        nrow=transformer.get_nrow();
        ncol=transformer.get_ncol();
        return;
    }

    delayed_reader(const delayed_reader& other) : original(other.original), 
            seed_ptr(other.seed_ptr->clone()), transformer(other.transformer) {}

    delayed_reader& operator=(const delayed_reader& other) {
        original=other.original;
        seed_ptr=other.seed_ptr->clone();
        transformer=other.transformer;
        return *this;
    }

    ~delayed_reader() = default;
    delayed_reader(delayed_reader&&) = default;
    delayed_reader& operator=(delayed_reader&&) = default;
    
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
    template<class Iter>
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    template<class Iter>
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    // Miscellaneous.
    Rcpp::RObject yield() const { return original; }

    static std::string get_class() { return "DelayedMatrix"; }

    static std::string get_package() { return "DelayedArray"; }
private:
    Rcpp::RObject original;
    std::unique_ptr<base_mat> seed_ptr;
    delayed_coord_transformer<T, V> transformer;

    // Specialized function for each realized matrix type.
    static std::unique_ptr<base_mat> generate_seed(Rcpp::RObject);
};

/******************************************************************
 * Implementing methods for the 'delayed_coord_transformer' class *
 ******************************************************************/

template<typename T, class V>
template<class M>
delayed_coord_transformer<T, V>::delayed_coord_transformer(M mat) : delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()) {}

template<typename T, class V>
template<class M>
delayed_coord_transformer<T, V>::delayed_coord_transformer(const Rcpp::List& net_subset, const Rcpp::LogicalVector& net_trans, M mat) :
        delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()), tmp(std::max(delayed_nrow, delayed_ncol)) { 
   
    const size_t original_nrow(mat->get_nrow()), original_ncol(mat->get_ncol()); 

    if (net_subset.size()!=2) {
        throw std::runtime_error("subsetting list should be of length 2");
    }

    // Checking indices for rows.
    obtain_indices(net_subset[0], original_nrow, byrow, delayed_nrow, row_index);

    // Checking indices for columns.
    obtain_indices(net_subset[1], original_ncol, bycol, delayed_ncol, col_index);

    // Checking transposition.
    if (net_trans.size()!=1) {
        throw std::runtime_error("transposition specifier should be of length 1");
    }
    transposed=net_trans[0];
    if (transposed) { // As the row/column indices refer to the matrix BEFORE transposition.
        std::swap(delayed_nrow, delayed_ncol);
    }

    return;
}

template<typename T, class V>
void delayed_coord_transformer<T, V>::obtain_indices(const Rcpp::RObject& subset_in, size_t original_dim,
        bool& affected, size_t& delayed_dim, std::vector<size_t>& subset_out) {
    // This function simply converts the row or column subset indices to 0-indexed form,
    // setting affected=false if there are no subset indices or if the subset indices are 1:original_dim.
    // Note that delayed_dim is also reset to the length fo the subset index vector.

    affected=!subset_in.isNULL();
    if (!affected){ 
        return;
    }

    // Coercing the subset indices to zero-indexed size_t's.
    if (subset_in.sexp_type()!=INTSXP) {
        throw std::runtime_error("index vector should be integer");
    }

    Rcpp::IntegerVector idx(subset_in);
    delayed_dim=idx.size();
    subset_out.reserve(delayed_dim);

    for (const auto i : idx) {
        if (i < 1 || static_cast<size_t>(i) > original_dim) {
            throw std::runtime_error("delayed subset indices are out of range");
        }
        subset_out.push_back(i-1);
    }

    // If the indices are all consecutive from 0 to N-1, we turn off 'affected'. 
    if (delayed_dim 
            && delayed_dim==original_dim
            && subset_out.front()==0 
            && subset_out.back()+1==delayed_dim) {

        size_t count=0;
        affected=false;
        for (auto i : subset_out) {
            if (i!=count) {
                affected=true;
                break;
            }
            ++count;
        }
    }

    return;
}

/*** Basic getter methods ***/

template<typename T, class V>
size_t delayed_coord_transformer<T, V>::get_nrow() const{ 
    return delayed_nrow;
}

template<typename T, class V>
size_t delayed_coord_transformer<T, V>::get_ncol() const{ 
    return delayed_ncol;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::get_row(M mat, size_t r, Iter out, size_t first, size_t last) {
    if (transposed) {
        dim_checker::check_dimension(r, get_nrow(), "row");
        dim_checker::check_subset(first, last, get_ncol(), "column");

        if (bycol) {
            r=col_index[r];
        }

        // Column extraction, first/last refer to rows.
        if (byrow) {
            reallocate_col(mat, r, first, last, out);
        } else {
            mat->get_col(r, out, first, last);
        }
    } else {
        if (byrow) {
            dim_checker::check_dimension(r, get_nrow(), "row");
            r=row_index[r];
        }

        // Row extraction, first/last refer to columns.
        if (bycol) {
            dim_checker::check_subset(first, last, get_ncol(), "column");
            reallocate_row(mat, r, first, last, out);
        } else {
            mat->get_row(r, out, first, last);
        }
    }
    return;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::get_col(M mat, size_t c, Iter out, size_t first, size_t last) {
    if (transposed) {
        dim_checker::check_dimension(c, get_ncol(), "column");
        dim_checker::check_subset(first, last, get_nrow(), "row");

        if (byrow) {
            c=row_index[c];
        }

        // Row extraction, first/last refer to columns.
        if (bycol) {
            reallocate_row(mat, c, first, last, out);
        } else {
            mat->get_row(c, out, first, last);
        }
    } else {
        if (bycol) {
            dim_checker::check_dimension(c, get_ncol(), "column");
            c=col_index[c];
        }

        // Column extraction, first/last refer to rows.
        if (byrow) {
            dim_checker::check_subset(first, last, get_nrow(), "row");
            reallocate_col(mat, c, first, last, out);
        } else {
            mat->get_col(c, out, first, last);
        }
    }
    return;
}

template<typename T, class V>
template<class M>
T delayed_coord_transformer<T, V>::get(M mat, size_t r, size_t c) {
    if (transposed) {
        dim_checker::check_dimension(r, get_nrow(), "row");
        dim_checker::check_dimension(c, get_ncol(), "column");
        if (bycol) {
            r=col_index[r];
        }
        if (byrow) {
            c=row_index[c];
        }
        return mat->get(c, r);
    } else {
        if (byrow) {
            dim_checker::check_dimension(r, get_nrow(), "row");
            r=row_index[r];
        }
        if (bycol) {
            dim_checker::check_dimension(c, get_ncol(), "column");
            c=col_index[c];
        }
        return mat->get(r, c);
    }
}

/*** Internal methods to handle reallocation after subsetting. ***/

template<typename T, class V>
void delayed_coord_transformer<T, V>::prepare_reallocation(size_t first, size_t last, 
        size_t& old_first, size_t& old_last, size_t& min_index, size_t& max_index, 
        const std::vector<size_t>& indices, const char* msg) {

    if (old_first!=first || old_last!=last) {
        old_first=first;
        old_last=last;
        if (first!=last) {
            min_index=*std::min_element(indices.begin()+first, indices.begin()+last);
            max_index=*std::max_element(indices.begin()+first, indices.begin()+last)+1;
        } else {
            // Avoid problems with max/min of zero-length vectors.
            min_index=0;
            max_index=0;
        }
    }

    return;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::reallocate_row(M mat, size_t r, size_t first, size_t last, Iter out) {
    // Yes, the use of *_col_* variables for row reallocation is intentional!
    // This is because we're talking about the columns in the extracted row.
    prepare_reallocation(first, last, old_col_first, old_col_last, 
            min_col_index, max_col_index, col_index, "column");

    V& holding=tmp.vec;
    mat->get_row(r, holding.begin(), min_col_index, max_col_index);
    auto cIt=col_index.begin()+first, end=col_index.begin()+last;
    while (cIt!=end) {
        (*out)=holding[*cIt - min_col_index];
        ++out;
        ++cIt;
    }
    return;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::reallocate_col(M mat, size_t c, size_t first, size_t last, Iter out) {
    // Yes, the use of *_row_* variables for column reallocation is intentional!
    // This is because we're talking about the rows in the extracted column.
    prepare_reallocation(first, last, old_row_first, old_row_last, 
            min_row_index, max_row_index, row_index, "row");

    V& holding=tmp.vec;
    mat->get_col(c, holding.begin(), min_row_index, max_row_index);
    auto rIt=row_index.begin()+first, end=row_index.begin()+last;
    while (rIt!=end) {
        (*out)=holding[*rIt - min_row_index];
        ++out;
        ++rIt;
    }
    return;
}

/*******************************************************
 * Implementing methods for the 'delayed_reader' class *
 *******************************************************/

/*** Constructor definitions ***/

/*** Basic getter methods ***/

template<typename T, class V, class base_mat>
template<class Iter>
void delayed_reader<T, V, base_mat>::get_col(size_t c, Iter out, size_t first, size_t last) {
    transformer.get_col(seed_ptr.get(), c, out, first, last);
    return;
}

template<typename T, class V, class base_mat>
template<class Iter>
void delayed_reader<T, V, base_mat>::get_row(size_t r, Iter out, size_t first, size_t last) {
    transformer.get_row(seed_ptr.get(), r, out, first, last);
    return;
}

template<typename T, class V, class base_mat>
T delayed_reader<T, V, base_mat>::get(size_t r, size_t c) {
    return transformer.get(seed_ptr.get(), r, c);
}

/*** Multi getter methods ***/

template<typename T, class V, class base_mat>
template<class Iter>
void delayed_reader<T, V, base_mat>::get_rows(Rcpp::IntegerVector::iterator rIt, size_t n, Iter out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);

    // No easy way to save memory or time, as we'd have to do a transposition if we use transformer.get_row().
    Rcpp::Environment beachenv(Rcpp::Environment::namespace_env("beachmat"));
    Rcpp::Function indexed_realizer(beachenv["realizeByIndexRange"]);

    Rcpp::IntegerVector cur_indices(rIt, rIt+n);
    for (auto& i : cur_indices) { ++i; }
    V tmp_store=indexed_realizer(original, cur_indices, Rcpp::IntegerVector::create(first, last-first));
    std::copy(tmp_store.begin(), tmp_store.end(), out);
    return;
}

template<typename T, class V, class base_mat>
template<class Iter>
void delayed_reader<T, V, base_mat>::get_cols(Rcpp::IntegerVector::iterator cIt, size_t n, Iter out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);

    if (seed_ptr->get_class()!="") { 
        // Not unknown, so there are probably fast column access methods available.
        const size_t nrows=last - first;
        for (size_t i=0; i<n; ++i, ++cIt, out+=nrows) {
            transformer.get_col(seed_ptr.get(), *cIt, out, first, last); 
        }
    } else {
        // Unknown matrices use block realization for speed (single block realization).
        Rcpp::Environment beachenv(Rcpp::Environment::namespace_env("beachmat"));
        Rcpp::Function indexed_realizer(beachenv["realizeByRangeIndex"]);

        Rcpp::IntegerVector cur_indices(cIt, cIt+n);
        for (auto& i : cur_indices) { ++i; }
        V tmp_store=indexed_realizer(original, Rcpp::IntegerVector::create(first, last-first), cur_indices);
        std::copy(tmp_store.begin(), tmp_store.end(), out);
    }
    return;
}

}

#endif
