#ifndef RATICATE_DELAYEDABIND_HPP
#define RATICATE_DELAYEDABIND_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <vector>
#include <memory>

namespace raticate {

template<typename Data = double, typename Index = int>
Parsed<Data, Index> parse_DelayedAbind(Rcpp::RObject seed) {
    Rcpp::List seeds(seed.slot("seeds"));

    Rcpp::List contents(seeds.size());
    std::vector<std::shared_ptr<tatami::Matrix<Data, Index> > > matrices(seeds.size());

    bool all_okay = true;
    for (size_t s = 0; s < seeds.size(); ++s) {
        auto sparsed = parse<Data, Index>(seeds[s]);
        if (sparsed.matrix == nullptr) {
            all_okay = false;
            break;
        }
        matrices[s] = sparsed.matrix;
        contents[s] = sparsed.contents;
    }

    Parsed<Data, Index> output;
    if (all_okay) {
        Rcpp::IntegerVector along(seed.slot("along"));
        if (along.size() != 1) {
            throw std::runtime_error("'along' should be an integer scalar");
        }

        if (along[0] == 1) {
            output.matrix = tatami::make_DelayedBind<0>(std::move(matrices));
        } else {
            output.matrix = tatami::make_DelayedBind<1>(std::move(matrices));
        }
        output.contents = contents;
    }

    return output;
}

}

#endif
