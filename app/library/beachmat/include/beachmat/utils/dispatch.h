#ifndef BEACHMAT_DISPATCH_H
#define BEACHMAT_DISPATCH_H

#include "Rcpp.h"

namespace beachmat {

template<class M>
std::unique_ptr<M> create_matrix(const Rcpp::RObject&);

template<class O>
std::unique_ptr<O> create_output(int, int, const output_param&);

};

#endif
