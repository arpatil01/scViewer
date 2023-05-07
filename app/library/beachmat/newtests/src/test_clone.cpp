#include "beachmat3/beachmat.h"
#include <algorithm>

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector test_clone(Rcpp::RObject mat) {
    auto thing = beachmat::read_lin_block(mat);
    auto thing2 = thing->clone();

    double value = 0;
    std::vector<double> workspace(thing->get_nrow());
    for (size_t i = 0; i < thing->get_ncol(); ++i) {
        if (i % 2 == 0) {
            auto out = thing->get_col(i, workspace.data());
            value += std::accumulate(out, out + thing->get_nrow(), 0.0);
        } else {
            auto out = thing2->get_col(i, workspace.data());
            value += std::accumulate(out, out + thing2->get_nrow(), 0.0);
        }
    }

    return Rcpp::NumericVector::create(value);
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector test_clone_sparse(Rcpp::RObject mat) {
    auto thing = beachmat::read_lin_sparse_block(mat);
    auto thing2 = thing->clone();

    double value = 0;
    std::vector<double> workspace_x(thing->get_nrow());
    std::vector<int> workspace_i(thing->get_nrow());

    for (size_t i = 0; i < thing->get_ncol(); ++i) {
        if (i % 2 == 0) {
            auto out = thing->get_col(i, workspace_x.data(), workspace_i.data());
            value += std::accumulate(out.x, out.x + out.n, 0.0);
        } else {
            auto out = thing2->get_col(i, workspace_x.data(), workspace_i.data());
            value += std::accumulate(out.x, out.x + out.n, 0.0);
        }
    }

    return Rcpp::NumericVector::create(value);
}
