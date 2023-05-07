#ifndef BEACHMAT_CSPARSE_READER_H
#define BEACHMAT_CSPARSE_READER_H

#include "Rcpp.h"

#include "../utils/utils.h"
#include "../utils/dim_checker.h"
#include "../utils/raw_structure.h"

#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace beachmat {

/*** Class definition ***/

template<typename T, class V>
class Csparse_reader : public dim_checker {
public:    
    Csparse_reader(const Rcpp::RObject&);
    ~Csparse_reader() = default;
    Csparse_reader(const Csparse_reader&) = default;
    Csparse_reader& operator=(const Csparse_reader&) = default;
    Csparse_reader(Csparse_reader&&) = default;
    Csparse_reader& operator=(Csparse_reader&&) = default;

    // Basic getters.
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
        in.n=get_const_col_nonzero(c, in.structure_start, in.values_start, first, last);
        return;
    }

    void get_row_raw(size_t r, raw_structure<V>& in, size_t first, size_t last) {
        check_rowargs(r, first, last);
        return;
    }

    static std::string col_raw_type () { return "sparse"; }

    static std::string row_raw_type () { return "none"; }

    // Multi getters.
    template <class Iter>
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);
    
    template <class Iter>
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    // Miscellaneous.
    Rcpp::RObject yield () const { return original; }

    static std::string get_class(); // specialized function for each realization.

    static std::string get_package() { return "Matrix"; }
protected:
    Rcpp::RObject original;
    Rcpp::IntegerVector i, p;
    V x;

    size_t currow, curstart, curend;
    std::vector<int> indices; // Left as 'int' to simplify comparisons with 'i' and 'p'.
    void update_indices(size_t, size_t, size_t);

    static T get_empty(); // Specialized function for each realization (easy to extend for non-int/double).

    size_t get_const_col_nonzero(size_t c, Rcpp::IntegerVector::iterator& index, typename V::iterator& val, size_t first, size_t last) {
        check_colargs(c, first, last);

        const int& pstart=p[c]; 
        index=i.begin()+pstart;
        auto endex=i.begin()+p[c+1]; 
        val=x.begin()+pstart;

        if (first) { // Jumping ahead if non-zero.
            auto new_index=std::lower_bound(index, endex, first);
            val+=(new_index-index);
            index=new_index;
        } 
        if (last!=(this->nrow)) { // Jumping to last element.
            endex=std::lower_bound(index, endex, last);
        }

        return endex-index;
    }
};

/*** Constructor definition ***/

template <typename T, class V>
Csparse_reader<T, V>::Csparse_reader(const Rcpp::RObject& incoming) : original(incoming), currow(0), curstart(0), curend(this->ncol) {
    auto classinfo=get_class_package(incoming);
    std::string ctype=classinfo.first;
    if (ctype!=get_class() || classinfo.second!=get_package()) {
        throw std::runtime_error(std::string("input should be a ") + ctype + " object");
    }

    this->fill_dims(get_safe_slot(incoming, "Dim"));
    const size_t& NC=this->ncol;
    const size_t& NR=this->nrow;

    Rcpp::RObject temp_i=get_safe_slot(incoming, "i");
    if (temp_i.sexp_type()!=INTSXP) { 
        throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should be integer"); 
    }
    i=temp_i;

    Rcpp::RObject temp_p=get_safe_slot(incoming, "p");
    if (temp_p.sexp_type()!=INTSXP) { 
        throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be integer");
    }
    p=temp_p;

    Rcpp::RObject temp_x=get_safe_slot(incoming, "x");
    if (temp_x.sexp_type()!=x.sexp_type()) { 
        throw std::runtime_error(std::string("'x' slot in a ") + ctype + " object should be " + translate_type(x.sexp_type()));
    }
    x=temp_x;

    if (x.size()!=i.size()) { 
        throw std::runtime_error(std::string("'x' and 'i' slots in a ") + ctype + " object should have the same length"); 
    }
    if (NC+1!=static_cast<size_t>(p.size())) { 
        throw std::runtime_error(std::string("length of 'p' slot in a ") + ctype + " object should be equal to 'ncol+1'"); 
    }
    if (p[0]!=0) { 
        throw std::runtime_error(std::string("first element of 'p' in a ") + ctype + " object should be 0"); 
    }
    if (p[NC]!=x.size()) { 
        throw std::runtime_error(std::string("last element of 'p' in a ") + ctype + " object should be 'length(x)'"); 
    }

    // Checking all the indices.
    auto pIt=p.begin();
    for (size_t px=0; px<NC; ++px) {
        const int& current=*pIt;
        if (current < 0) { 
            throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should contain non-negative values"); 
        }
        if (current > *(++pIt)) { 
            throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be sorted"); 
        }
    }

    pIt=p.begin();
    for (size_t px=0; px<NC; ++px) {
        int left=*pIt; // Integers as that's R's storage type. 
        int right=*(++pIt)-1; // Not checking the last element, as this is the start of the next column.
        auto iIt=i.begin()+left;

        for (int ix=left; ix<right; ++ix) {
            const int& current=*iIt;
            if (current > *(++iIt)) {
                throw std::runtime_error(std::string("'i' in each column of a ") + ctype + " object should be sorted");
            }
        }
    }

    for (auto iIt=i.begin(); iIt!=i.end(); ++iIt) {
        const int& curi=*iIt;
        if (curi<0 || static_cast<size_t>(curi)>=NR) {
            throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should contain elements in [0, nrow)");
        }
    }

    return;
}

/*** Basic getter functions ***/

template <typename T, class V>
T Csparse_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    auto iend=i.begin() + p[c+1];
    auto loc=std::lower_bound(i.begin() + p[c], iend, r);
    if (loc!=iend && static_cast<size_t>(*loc)==r) { 
        return x[loc - i.begin()];
    } else {
        return get_empty();
    }
}

template <typename T, class V>
void Csparse_reader<T, V>::update_indices(size_t r, size_t first, size_t last) {
    /* Initializing the indices upon the first request, assuming currow=0 based on initialization above.
     * This avoids using up space for the indices if we never do row access.
     */
    if (indices.size()!=this->ncol) {
        indices=std::vector<int>(p.begin(), p.begin()+this->ncol);
    }

    /* If left/right slice are not equal to what is stored, we reset the indices,
     * so that the code below will know to recompute them. It's too much effort
     * to try to figure out exactly which columns need recomputing; just do them all.
     */
    if (first!=curstart || last!=curend) {
        curstart=first;
        curend=last;
        Rcpp::IntegerVector::iterator pIt=p.begin()+first;
        for (size_t px=first; px<last; ++px, ++pIt) {
            indices[px]=*pIt; 
        }
        currow=0;
    }

    /* entry of 'indices' for each column should contain the index of the first
     * element with row number not less than 'r'. If no such element exists, it
     * will contain the index of the first element of the next column.
     */
    if (r==currow) { 
        return; 
    } 

    Rcpp::IntegerVector::iterator pIt=p.begin()+first;
    if (r==currow+1) {
        ++pIt; // points to the first-past-the-end element, at any given 'c'.
        for (size_t c=first; c<last; ++c, ++pIt) {
            int& curdex=indices[c];
            if (curdex!=*pIt && static_cast<size_t>(i[curdex]) < r) { 
                ++curdex;
            }
        }
    } else if (r+1==currow) {
        for (size_t c=first; c<last; ++c, ++pIt) {
            int& curdex=indices[c];
            if (curdex!=*pIt && static_cast<size_t>(i[curdex-1]) >= r) { 
                --curdex;
            }
        }

    } else { 
        Rcpp::IntegerVector::iterator istart=i.begin(), loc;
        if (r > currow) {
            ++pIt; // points to the first-past-the-end element, at any given 'c'.
            for (size_t c=first; c<last; ++c, ++pIt) { 
                int& curdex=indices[c];
                loc=std::lower_bound(istart + curdex, istart + *pIt, r);
                curdex=loc - istart;
            }
        } else { 
            for (size_t c=first; c<last; ++c, ++pIt) {
                int& curdex=indices[c];
                loc=std::lower_bound(istart + *pIt, istart + curdex, r);
                curdex=loc - istart;
            }
        }
    }

    currow=r;
    return;
}

template <typename T, class V>
template <class Iter>
void Csparse_reader<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    update_indices(r, first, last);
    std::fill(out, out+last-first, get_empty());

    auto pIt=p.begin()+first+1; // Points to first-past-the-end for each 'c'.
    for (size_t c=first; c<last; ++c, ++pIt, ++out) { 
        const int& idex=indices[c];
        if (idex!=*pIt && static_cast<size_t>(i[idex])==r) { (*out)=x[idex]; }
    } 
    return;  
}

template <typename T, class V>
template <class Iter>
void Csparse_reader<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    const int& pstart=p[c]; 
    auto iIt=i.begin()+pstart, 
         eIt=i.begin()+p[c+1]; 
    auto xIt=x.begin()+pstart;

    if (first) { // Jumping ahead if non-zero.
        auto new_iIt=std::lower_bound(iIt, eIt, first);
        xIt+=(new_iIt-iIt);
        iIt=new_iIt;
    } 
    if (last!=(this->nrow)) { // Jumping to last element.
        eIt=std::lower_bound(iIt, eIt, last);
    }

    std::fill(out, out+last-first, get_empty());
    for (; iIt!=eIt; ++iIt, ++xIt) {
        *(out + (*iIt - int(first)))=*xIt;
    }
    return;
}

/*** Multi getter functions ***/

template<typename T, class V>
template<class Iter>
void Csparse_reader<T, V>::get_rows(Rcpp::IntegerVector::iterator rIt, size_t n, Iter out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);

    Rcpp::IntegerVector::iterator indices;
    typename V::iterator values;

    for (size_t c=first; c<last; ++c) {
		size_t nnzero=get_const_col_nonzero(c, indices, values, 0, this->nrow);
        auto endpoint=indices+nnzero;
        auto rIt_copy=rIt;

        for (size_t i=0; i<n; ++i, ++out, ++rIt_copy) {  
            if (indices==endpoint) {
                (*out)=0;
            } else if (*rIt_copy==*indices) {
                (*out)=*values;
                ++indices;
                ++values;
            } else if (*rIt_copy < *indices) {
                (*out)=0;
            } else {
                auto next=std::lower_bound(indices, endpoint, *rIt_copy);
                values+=(next - indices);
                indices=next;
                
                if (next!=endpoint && *next==*rIt_copy) {
                    (*out)=*values;
                    ++indices;
                    ++values;
                } else {
                    (*out)=0;
                }
            }
        }
    }
    return;
}

template<typename T, class V>
template<class Iter>
void Csparse_reader<T, V>::get_cols(Rcpp::IntegerVector::iterator cIt, size_t n, Iter out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);
    size_t nrows=last - first;
    for (size_t i=0; i<n; ++i, ++cIt, out+=nrows) {
        get_col(*cIt, out, first, last);
    }
    return;
}

}

#endif
