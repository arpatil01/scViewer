#include "beachmat3/beachmat.h"
#include <algorithm>
#include <vector>

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector test_promotion(Rcpp::RObject mat) {
    auto ptr = beachmat::read_lin_block(mat);
    double value = 0;

    if (ptr->is_sparse()) {
        auto sptr = beachmat::promote_to_sparse(ptr);
        std::vector<double> workspace_x(sptr->get_nrow());
        std::vector<int> workspace_i(sptr->get_nrow());

        for (size_t i = 0; i < sptr->get_ncol(); ++i) {
            auto out = sptr->get_col(i, workspace_x.data(), workspace_i.data());
            value += std::accumulate(out.x, out.x + out.n, 0.0);
        }
    } else {
        std::vector<double> workspace(ptr->get_nrow());
        for (size_t i = 0; i < ptr->get_ncol(); ++i) {
            auto out = ptr->get_col(i, workspace.data());
            value += std::accumulate(out, out + ptr->get_nrow(), 0.0);
        }

        ++value; // as a mark that we hit this code clause.
    }

    return Rcpp::NumericVector::create(value);
}
