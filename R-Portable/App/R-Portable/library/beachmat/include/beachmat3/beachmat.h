#ifndef BEACHMAT_BEACHMAT_H
#define BEACHMAT_BEACHMAT_H

/**
 * @mainpage 
 *
 * The **beachmat** R package provides a C++ API for extracting row/column vectors from dense/sparse blocks of a matrix, 
 * typically generated in each iteration of the **DelayedArray** block processing mechanism.
 * Readers are referred to http://bioconductor.org/packages/devel/bioc/html/beachmat.html
 * for a more user-friendly introduction into the purpose and usage of **beachmat**;
 * this site simply provides the reference documentation for the various C++ classes and functions.
 *
 * @authors Aaron Lun
 */

/**
 * @file beachmat.h
 *
 * The main header file that should be `include`d to use the **beachmat** API.
 * This provides functions to easily construct an instance of the appropriate subclass from an R object.
 */

#include "Rcpp.h"
#include "read_lin_block.h"
#include "as_gCMatrix.h"
#include <stdexcept>
#include <memory>

#endif
