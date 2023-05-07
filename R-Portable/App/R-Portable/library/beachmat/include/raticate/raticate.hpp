#ifndef RATICATE_HPP
#define RATICATE_HPP

#include "Rcpp.h"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include "CSparseMatrix.hpp"
#include "DelayedMatrix.hpp"
#include "DelayedSubset.hpp"
#include "DelayedAperm.hpp"
#include "DelayedAbind.hpp"
#include "UnknownMatrix.hpp"
#include "utils.hpp"
#include "parallelize.hpp"

/**
 * @file raticate.hpp
 *
 * @brief Parse `Rcpp::RObject`s into **tatami** matrices.
 */

namespace raticate {

/**
 * Parse `Rcpp::RObject`s into **tatami** matrices.
 * Natively supported matrix types are:
 *
 * - ordinary logical, numeric or integer matrices.
 * - `dgCMatrix` or `lgCMatrix` objects from the **Matrix** package.
 * - `SparseArraySeed` objects from the [**DelayedArray**](https://github.com/Bioconductor/DelayedArray) package.
 * - `DelayedMatrix` objects wrapping any of the above, or containing the following delayed operations:
 *    - Subsetting (as a `DelayedSubset` instance)
 *    - Modification of dimnames (as a `DelayedSetDimnames` instance)
 *    - Transposition (as a `DelayedAperm` instance)
 *    - Combining (as a `DelayedAbind` instance)
 * 
 * For all other objects, if `allow_unknown = true`, we create an instance of an "unknown matrix fallback" subclass of a `tatami::Matrix`.
 * This calls `DelayedArray::extract_array()` in R to extract an appropriate slice of the matrix.
 * Of course, this is quite a bit slower than the native representations.
 * Also see [here](../../docs/parallel.md) for safe parallelization of **tatami** calls that might operate on an unknown matrix.
 *
 * If `allow_unknown = false`, any object not listed above will result in a `nullptr`.
 * This should be handled by the caller. 
 * 
 * @tparam Data Numeric data type for the **tatami** interface, typically `double`.
 * @tparam Index Integer index type for the **tatami** interface, typically `int`.
 * 
 * @param x An R object representing a supported matrix type.
 * @param allow_unknown Whether to allow the creation of a fallback `tatami::Matrix` for unknown matrix-like objects. 
 *
 * @return A `Parsed` object containing a pointer to a parsed `tatami::Matrix` (or `null`, if parsing was not successful and `allow_unknown = false`).
 */
template<typename Data, typename Index>
Parsed<Data, Index> parse(Rcpp::RObject x, bool allow_unknown /* = false, in the declaration at utils.hpp */) {
    Parsed<Data, Index> output;

    if (x.isS4()) {
        std::string ctype = get_class_name(x);

        if (ctype == "SparseArraySeed") {
            output = parse_SparseArraySeed<Data, Index>(x);
        } else if (ctype == "dgCMatrix") {
            output = parse_dgCMatrix<Data, Index>(x);
        } else if (ctype == "lgCMatrix") {
            output = parse_lgCMatrix<Data, Index>(x);

        } else if (ctype == "DelayedMatrix") {
            output = parse_DelayedMatrix<Data, Index>(x);
        } else if (ctype == "DelayedSetDimnames") {
            output = parse_DelayedMatrix<Data, Index>(x); // just forward onto the seed.
        } else if (ctype == "DelayedSubset") {
            output = parse_DelayedSubset<Data, Index>(x);
        } else if (ctype == "DelayedAperm") {
            output = parse_DelayedAperm<Data, Index>(x);
        } else if (ctype == "DelayedAbind") {
            output = parse_DelayedAbind<Data, Index>(x);
        }

    } else if (x.hasAttribute("dim")) {
        output = parse_simple_matrix<Data, Index>(x);
    }

    /**
     * As a general rule, any internal calls to parse() in DelayedArray parsers
     * should keep the default of allow_unknown = false. This avoids partial
     * parsing that ends up having to fall back to R anyway. By keeping as much
     * as we can in R, we might be able to get more effective extraction via
     * block processing, which sees multiple indices at once for better
     * extraction (tatami only sees one row/column at a time).
     */

    if (output.matrix == nullptr && allow_unknown) {
        // No need to set contents here, as the matrix itself holds the Rcpp::RObject.
        output.matrix.reset(new UnknownMatrix<Data, Index>(x));
    }

    return output;
}

}

#endif
