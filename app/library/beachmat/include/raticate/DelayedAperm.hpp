#ifndef RATICATE_DELAYEDAPERM_HPP
#define RATICATE_DELAYEDAPERM_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <vector>

namespace raticate {

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_DelayedAperm(Rcpp::RObject seed) {
    auto sparsed = parse<Data, Index>(seed.slot("seed"));

    if (sparsed.matrix != nullptr) {
        Rcpp::IntegerVector perm(seed.slot("perm"));
        if (perm.size() != 2) {
            throw std::runtime_error("'perm' slot should be an integer vector of length 2");
        }

        // Seeing if we actually need to permute it.
        if (perm[0] == 2 && perm[1] == 1) {
            sparsed.matrix = tatami::make_DelayedTranspose(sparsed.matrix);
        }
    }

    return sparsed;
}

}

#endif
