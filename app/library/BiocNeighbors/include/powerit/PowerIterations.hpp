#ifndef POWERIT_POWERIT_HPP
#define POWERIT_POWERIT_HPP

#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include "aarand/aarand.hpp"

/**
 * @file powerit.hpp
 *
 * @brief Implements some power iterations.
 */

namespace powerit {

/**
 * @brief Perform power iterations to find the first eigenvalue/vector.
 *
 * Pretty much as it says. 
 * Given a diagonalizable matrix, this will perform power iterations to find the first eigenvalue and vector.
 * For example, we can apply this to a covariance matrix to obtain the first principal component.
 */
class PowerIterations {
public:
    /**
     * @brief Default parameter settings.
     */
    struct Defaults {
        /**
         * See `set_iterations()` for details.
         */
        static constexpr int iterations = 500;

        /**
         * See `set_tolerance()` for details.
         */
        static constexpr double tolerance = 0.000001;
    };

    /**
     * @param i Maximum number of iterations to perform.
     * Note that the algorithm may converge before reaching this limit.
     *
     * @return A reference to this `PowerIterations` object.
     */
    PowerIterations& set_iterations(int i = Defaults::iterations) {
        iterations = i;
        return *this;
    }

    /**
     * @param t Tolerance used to determine convergence.
     * The error is defined as the L2 norm of the difference between eigenvectors at successive iterations;
     * if this drops below `t`, we assume that the algorithm has converged.
     *
     * @return A reference to this `PowerIterations` object.
     */
    PowerIterations& set_tolerance(double t = Defaults::tolerance) {
        tolerance = t;
        return *this;
    }

private:
    int iterations = Defaults::iterations;
    double tolerance = Defaults::tolerance;

    template<typename Data>
    static Data normalize(int ndim, Data* x) {
        Data ss = 0;
        for (int d = 0; d < ndim; ++d) {
            ss += x[d] * x[d];
        }

        if (ss) {
            ss = std::sqrt(ss);
            for (int d = 0; d < ndim; ++d) {
                x[d] /= ss;
            }
        }
        return ss;
    }

public:
    /**
     * @tparam Data Floating-point type for the data.
     * @tparam Engine Any C++11-compliant random number generator class.
     *
     * @param order Order of the square matrix.
     * @param[in] x Pointer to an array containing an `order`-by-`order` diagonalizable matrix.
     * This is assumed to be column-major.
     * @param[out] output Pointer to an array of length `order`,
     * to store the output eigenvector.
     * @param eng Instance of the random number generator.
     *
     * @return Pair containing the eigenvalue (first) and the number of iterations required for convergence (second).
     * The latter is set to -1 if convergence did not occur before the maximum number of iterations specified in `set_iterations()`.
     */
    template<class Data, class Engine>
    std::pair<Data, int> run(size_t order, const Data* x, Data* output, Engine& eng) {
        // Defining a random starting vector.
        while (1) {
            for (int d = 0; d < order - 1; d += 2) {
                auto sampled = aarand::standard_normal<Data>(eng);
                output[d] = sampled.first;
                output[d + 1] = sampled.second;
            }
            if (order % 2) {
                output[order - 1] = aarand::standard_normal<Data>(eng).first;
            }
            if (normalize(order, output)) {
                break;
            }
        }

        std::pair<Data, int> stats;
        auto& l2 = stats.first;
        stats.second = -1;
        std::vector<Data> buffer(order);

        // Applying power iterations.
        for (int i = 0; i < iterations; ++i) {
            for (size_t j = 0; j < order; ++j) {
                // As the matrix is symmetric, we can use inner_product.
                // This technically computes the transpose of the matrix
                // with the vector, but it's all the same, so whatever.
                buffer[j] = std::inner_product(output, output + order, x + j * order, static_cast<Data>(0.0));
            }

            // Normalizing the vector.
            l2 = normalize(order, buffer.data());

            // We want to know if SIGMA * output = lambda * output, i.e., l2 * buffer = lambda * output.
            // If we use l2 as a working estimate for lambda, we're basically just testing the difference
            // between buffer and output. We compute the error and compare this to the tolerance. No need
            // to normalize as we've already implicitly done that.
            Data err = 0;
            for (size_t d = 0; d < order; ++d) {
                Data diff = buffer[d] - output[d];
                err += diff * diff;
            }
            if (std::sqrt(err) < tolerance) {
                stats.second = i + 1;
                break;
            }

            std::copy(buffer.begin(), buffer.end(), output);
        }

        return stats;
    } 
};

}

#endif
