#include "beachmat3/beachmat.h"
#include <algorithm>

template <class M, typename T = typename M::stored_type>
Rcpp::RObject get_row_slice0(Rcpp::RObject mat, Rcpp::IntegerVector order,
    Rcpp::IntegerVector starts, Rcpp::IntegerVector ends)
{
    auto ptr = beachmat::read_lin_block(mat);
    std::vector<T> tmp(ptr->get_ncol());
    M output(ptr->get_nrow(), ptr->get_ncol());

    for (auto o : order) {
        int curstart = starts[o];
        int curend = ends[o];
        auto vec = ptr->get_row(o, tmp.data(), curstart, curend);
        auto curout = output.row(o);
        std::copy(vec, vec + curend - curstart, curout.begin() + curstart);
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject get_row_slice(Rcpp::RObject mat, Rcpp::IntegerVector order,
    Rcpp::IntegerVector starts, Rcpp::IntegerVector ends, int mode)
{
    if (mode==0) {
        return get_row_slice0<Rcpp::LogicalMatrix>(mat, order, starts, ends);
    } else if (mode==1) {
        return get_row_slice0<Rcpp::IntegerMatrix>(mat, order, starts, ends);
    } else {
        return get_row_slice0<Rcpp::NumericMatrix>(mat, order, starts, ends);
    }
}

template <class M, typename T = typename M::stored_type>
Rcpp::RObject get_row0(Rcpp::RObject mat, Rcpp::IntegerVector order) {
    auto ptr = beachmat::read_lin_block(mat);
    std::vector<T> tmp(ptr->get_ncol());
    M output(ptr->get_nrow(), ptr->get_ncol());

    for (auto o : order) {
        auto vec = ptr->get_row(o, tmp.data());
        auto curout = output.row(o);
        std::copy(vec, vec + ptr->get_ncol(), curout.begin());
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject get_row(Rcpp::RObject mat, Rcpp::IntegerVector order, int mode)
{
    if (mode==0) {
        return get_row0<Rcpp::LogicalMatrix>(mat, order);
    } else if (mode==1) {
        return get_row0<Rcpp::IntegerMatrix>(mat, order);
    } else {
        return get_row0<Rcpp::NumericMatrix>(mat, order);
    }
}
