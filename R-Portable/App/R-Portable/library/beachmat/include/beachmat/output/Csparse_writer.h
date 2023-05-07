#ifndef BEACHMAT_CSPARSE_WRITER_H
#define BEACHMAT_CSPARSE_WRITER_H

#include "Rcpp.h"

#include "../utils/utils.h"
#include "../utils/dim_checker.h"

#include <utility>
#include <vector>
#include <deque>
#include <algorithm>

namespace beachmat { 

/*** Class definition ***/

template<typename T, class V>
class Csparse_writer : public dim_checker {
public:
    Csparse_writer(size_t, size_t);
    ~Csparse_writer() = default;
    Csparse_writer(const Csparse_writer&) = default;
    Csparse_writer& operator=(const Csparse_writer&) = default;
    Csparse_writer(Csparse_writer&&) = default;
    Csparse_writer& operator=(Csparse_writer&&) = default;

    // Setters:
    template <class Iter>
    void set_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void set_col(size_t, Iter, size_t, size_t);

    void set(size_t, size_t, T);
    
    template <class Iter>
    void set_col_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Iter);

    template <class Iter>
    void set_row_indexed(size_t, size_t, Rcpp::IntegerVector::iterator, Iter);

    // Getters:
    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_row(size_t, Iter, size_t, size_t);

    T get(size_t, size_t);

    // Other:
    Rcpp::RObject yield();

    static std::string get_class();

    static std::string get_package() { return "Matrix"; }
private:
    typedef std::pair<size_t, T> data_pair;
    std::vector<std::deque<data_pair> > data;

    // What is an empty value?
    static T get_empty();

    // Only comparing the first value.
    static bool only_first_less(const data_pair& lhs, const data_pair& rhs) { return lhs.first < rhs.first; }

    template <class Iter>
    static Iter find_matching_row(Iter, Iter, const data_pair&);

    // General column insertions.
    static void insert_into_column(std::deque<data_pair>&, size_t, T);
};

/*** Constructor definition ***/

template<typename T, class V>
Csparse_writer<T, V>::Csparse_writer(size_t nr, size_t nc) : dim_checker(nr, nc), data(nc) {}

/*** Setter methods ***/

template<typename T, class V>
template<class Iter>
void Csparse_writer<T, V>::set_col(size_t c, Iter in, size_t first, size_t last) {
    check_colargs(c, first, last);
    std::deque<data_pair>& current=data[c];
    std::deque<data_pair> new_set;

    // Filling in all elements before start.
    auto cIt=current.begin(); 
    while (cIt!=current.end() && cIt->first < first) {
        new_set.push_back(*cIt);
        ++cIt;
    }
   
    // Filling in all non-empty elements. 
    for (size_t index=first; index<last; ++index, ++in) {
        if ((*in)!=get_empty()) { 
            new_set.push_back(data_pair(index, *in));
        }
    } 

    // Jumping to the end.
    while (cIt!=current.end() && cIt->first < last) {
        ++cIt;
    }

    // Filling in remaining elements.
    while (cIt!=current.end()) {
        new_set.push_back(*cIt);
        ++cIt;
    }

    current.swap(new_set);
    return;
}

template<typename T, class V>
template <class Iter>
Iter Csparse_writer<T, V>::find_matching_row(Iter begin, Iter end, const data_pair& incoming) {
    return std::lower_bound(begin, end, incoming, only_first_less);
}

template<typename T, class V>
void Csparse_writer<T, V>::insert_into_column(std::deque<data_pair>& column, size_t r, T val) {
    if (column.size()) {
        if (r < column.front().first) {
            column.push_front(data_pair(r, val));
        } else if (r==column.front().first) {
            column.front().second=val;
        } else if (r > column.back().first) {
            column.push_back(data_pair(r, val));
        } else if (r==column.back().first) {
            column.back().second=val;
        } else {
            data_pair incoming(r, val);
            auto insert_loc=find_matching_row(column.begin(), column.end(), incoming);
            if (insert_loc!=column.end() && insert_loc->first==r) { 
                insert_loc->second=val;
            } else {
                column.insert(insert_loc, incoming);
            }
        }
    } else {
        column.push_back(data_pair(r, val));
    }
    return;
}

template<typename T, class V>
template<class Iter>
void Csparse_writer<T, V>::set_row(size_t r, Iter in, size_t first, size_t last) {
    check_rowargs(r, first, last);
    for (size_t c=first; c<last; ++c, ++in) {
        if ((*in)==get_empty()) { continue; }
        insert_into_column(data[c], r, *in);
    }
    return;
}

template<typename T, class V>
void Csparse_writer<T, V>::set(size_t r, size_t c, T in) {
    check_oneargs(r, c);
    set_row(r, &in, c, c+1);
    return;
}

template<typename T, class V>
template <class Iter>
void Csparse_writer<T, V>::set_col_indexed(size_t c, size_t n, Rcpp::IntegerVector::iterator idx, Iter in) {
    check_colargs(c);
    std::deque<data_pair>& current=data[c];
    for (size_t i=0; i<n; ++i, ++idx, ++in) { 
        current.push_back(data_pair(*idx, *in));
    }

    std::stable_sort(current.begin(), current.end(), only_first_less);
    std::deque<data_pair> survivors;
    auto cIt=current.begin();
    while (cIt!=current.end()) { 
        auto first=cIt->first; 
        ++cIt;
        while (cIt!=current.end() && cIt->first==first) { ++cIt; }
        survivors.push_back(*(cIt-1));
    }

    current.swap(survivors);
    return;
}

template<typename T, class V>
template <class Iter>
void Csparse_writer<T, V>::set_row_indexed(size_t r, size_t n, Rcpp::IntegerVector::iterator idx, Iter in) {
    check_rowargs(r);
    for (size_t i=0; i<n; ++i, ++idx, ++in) { 
        insert_into_column(data[*idx], r, *in);
    }
    return;
}

/*** Getter methods ***/

template<typename T, class V>
template<class Iter>
void Csparse_writer<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    // It is not easy to use caching here, like Csparse_matrix() does. This is
    // because cached indices can be invalidated upon set_row().
    check_rowargs(r, first, last);
    std::fill(out, out+last-first, get_empty());

    for (size_t col=first; col<last; ++col, ++out) {
        const std::deque<data_pair>& current=data[col];
        if (current.empty() || r>current.back().first || r<current.front().first) {
            continue; 
        }
        if (r==current.back().first) { 
            (*out)=current.back().second;
        } else if (r==current.front().first) {
            (*out)=current.front().second;
        } else {
            // (*out) is equivalent to get_empty(), due to fill; not that it matters, 
            // as find_matching_row() will only use the first value anyway.
            auto loc=find_matching_row(current.begin(), current.end(), data_pair(r, *out)); 
            if (loc!=current.end() && loc->first==r) { 
                (*out)=loc->second;
            }
        }
    }
    return;
}

template<typename T, class V>
template<class Iter>
void Csparse_writer<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    const std::deque<data_pair>& current=data[c];

    // Jumping forwards.
    auto cIt=current.begin();
    if (first) {
        cIt=find_matching_row(current.begin(), current.end(), data_pair(first, get_empty()));
    }
    
    std::fill(out, out+last-first, get_empty());
    while (cIt!=current.end() && cIt->first < last) { 
        *(out + (cIt->first - first)) = cIt->second;
        ++cIt;
    }
    return;
}

template<typename T, class V>
T Csparse_writer<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    const std::deque<data_pair>& current=data[c];
    auto cIt=find_matching_row(current.begin(), current.end(), data_pair(r, get_empty()));
    if (cIt!=current.end() && cIt->first==r) {
        return cIt->second;
    } else {
        return get_empty();
    }
}

/*** Output function ***/

template<typename T, class V>
Rcpp::RObject Csparse_writer<T, V>::yield() {
    std::string classname=get_class();
    Rcpp::S4 mat(classname);

    // Setting dimensions.
    if (!mat.hasSlot("Dim")) {
        throw std::runtime_error(std::string("missing 'Dim' slot in ") + classname + " object");
    }
    mat.slot("Dim") = Rcpp::IntegerVector::create(this->nrow, this->ncol);

    // Setting 'p'.
    if (!mat.hasSlot("p")) {
        throw std::runtime_error(std::string("missing 'p' slot in ") + classname + " object");
    }
    Rcpp::IntegerVector p(this->ncol+1, 0);
    auto pIt=p.begin()+1;
    size_t total_size=0;
    for (auto dIt=data.begin(); dIt!=data.end(); ++dIt, ++pIt) { 
        total_size+=dIt->size();
        (*pIt)=total_size;
    }
    mat.slot("p")=p;

    // Setting 'i' and 'x'.
    Rcpp::IntegerVector i(total_size);
    V x(total_size);
    if (!mat.hasSlot("i")) {
        throw std::runtime_error(std::string("missing 'i' slot in ") + classname + " object");
    }
    if (!mat.hasSlot("x")) {
        throw std::runtime_error(std::string("missing 'x' slot in ") + classname + " object");
    }
    auto xIt=x.begin();
    auto iIt=i.begin();
    for (size_t c=0; c<this->ncol; ++c) {
        auto current=data[c];
        for (auto cIt=current.begin(); cIt!=current.end(); ++cIt, ++xIt, ++iIt) {
            (*iIt)=cIt->first;
            (*xIt)=cIt->second;
        }
    }
    mat.slot("i")=i;
    mat.slot("x")=x;

    return SEXP(mat);
}

}

#endif
