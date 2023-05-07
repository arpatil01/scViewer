#ifndef RATICATE_DELAYEDMATRIX_HPP
#define RATICATE_DELAYEDMATRIX_HPP

#include "utils.hpp"

namespace raticate {

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_DelayedMatrix(Rcpp::RObject seed) {
    return parse<Data, Index>(seed.slot("seed"));
}

}

#endif
