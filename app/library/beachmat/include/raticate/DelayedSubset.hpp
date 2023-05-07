#ifndef RATICATE_DELAYEDSUBSET_HPP
#define RATICATE_DELAYEDSUBSET_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <vector>

namespace raticate {

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_DelayedSubset(Rcpp::RObject seed) {
    auto sparsed = parse<Data, Index>(seed.slot("seed"));
    auto mat = sparsed.matrix;

    if (mat != nullptr) {
        Rcpp::List index(seed.slot("index"));
        if (index.size() != 2) {
            throw std::runtime_error("'index' slot should be a list of length 2");
        }

        // Seeing if we need to subset by row.
        Rcpp::RObject byrow(index[0]);
        if (!byrow.isNULL()) {
            Rcpp::IntegerVector rindex(byrow);
            std::vector<int> copy(rindex.begin(), rindex.end());
            for (auto& c : copy) { --c; } // 1-based indexing.

            mat = tatami::make_DelayedSubset<0>(mat, std::move(copy));       
        }

        // Seeing if we need to subset by column.
        Rcpp::RObject bycol(index[1]);
        if (!bycol.isNULL()) {
            Rcpp::IntegerVector cindex(bycol);
            std::vector<int> copy(cindex.begin(), cindex.end());
            for (auto& c : copy) { --c; } // 1-based indexing.

            mat = tatami::make_DelayedSubset<1>(mat, std::move(copy));       
        }

        sparsed.matrix = mat;
    }

    return sparsed;
}

}

#endif
