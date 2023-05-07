#ifndef RATICATE_PARALLELIZE_HPP
#define RATICATE_PARALLELIZE_HPP

/**
 * @cond
 */
#ifdef RATICATE_PARALLELIZE_UNKNOWN
/**
 * @endcond
 */

#include "UnknownMatrix.hpp"

/**
 * @file parallelize.hpp
 *
 * @brief Safely parallelize for unknown matrix fallbacks.
 */

namespace raticate {

/**
 * @tparam Data Numeric data type for the **tatami** interface, typically `double`.
 * @tparam Index Integer index type for the **tatami** interface, typically `int`.
 *
 * @param njobs Number of jobs to be executed.
 * @param fun Function to run in each thread.
 * This is typically a lambda that should accept two arguments specifying the first and one-past-the-last job to be executed in a given thread.
 * @param nthreads Number of threads to parallelize over.
 *
 * The series of integers from 0 to `njobs - 1` is split into `nthreads` contiguous ranges.
 * Each range is used as input to `fun` within the corresponding thread.
 * It is assumed that the execution of any given job is independent of the next.
 */ 
template<typename Data, typename Index, class Function>
void parallelize(size_t njobs, Function fun, size_t nthreads) {
    parallel_coordinator().template run<Data, Index>(njobs, fun, nthreads);
}

}

/**
 * @cond
 */
#endif
/**
 * @endcond
 */

#endif
