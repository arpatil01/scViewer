#include "beachmat3/beachmat.h"
#include <algorithm>

template <class M, typename T = typename M::stored_type>
Rcpp::RObject get_column_slice0(Rcpp::RObject mat, Rcpp::IntegerVector order, 
    Rcpp::IntegerVector starts, Rcpp::IntegerVector ends) 
{
    auto ptr = beachmat::read_lin_block(mat);
    std::vector<T> tmp(ptr->get_nrow());
    M output(ptr->get_nrow(), ptr->get_ncol());

    for (auto o : order) {
        int curstart = starts[o];
        int curend = ends[o];
        auto vec = ptr->get_col(o, tmp.data(), curstart, curend);

        auto curout = output.column(o);
        std::copy(vec, vec + curend - curstart, curout.begin() + curstart);
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject get_column_slice(Rcpp::RObject mat, Rcpp::IntegerVector order, 
    Rcpp::IntegerVector starts, Rcpp::IntegerVector ends, int mode) 
{
    if (mode==0) {
        return get_column_slice0<Rcpp::LogicalMatrix>(mat, order, starts, ends);
    } else if (mode==1) {
        return get_column_slice0<Rcpp::IntegerMatrix>(mat, order, starts, ends);
    } else {
        return get_column_slice0<Rcpp::NumericMatrix>(mat, order, starts, ends);
    }
}

template <class M, typename T = typename M::stored_type>
Rcpp::RObject get_column0(Rcpp::RObject mat, Rcpp::IntegerVector order) {
    auto ptr = beachmat::read_lin_block(mat);
    std::vector<T> tmp(ptr->get_nrow());
    M output(ptr->get_nrow(), ptr->get_ncol());

    for (auto o : order) {
        auto vec = ptr->get_col(o, tmp.data());
        auto curout = output.column(o);
        std::copy(vec, vec + ptr->get_nrow(), curout.begin());
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject get_column(Rcpp::RObject mat, Rcpp::IntegerVector order, int mode)
{
    if (mode==0) {
        return get_column0<Rcpp::LogicalMatrix>(mat, order);
    } else if (mode==1) {
        return get_column0<Rcpp::IntegerMatrix>(mat, order);
    } else {
        return get_column0<Rcpp::NumericMatrix>(mat, order);
    }
}
