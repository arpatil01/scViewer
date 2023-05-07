#ifndef BEACHMAT_UTILS_H
#define BEACHMAT_UTILS_H

/**
 * @file utils.h
 *
 * Internal utilities for use in the various **beachmat** classes.
 */

#include "Rcpp.h"
#include <string>
#include <utility>
#include <stdexcept>

namespace beachmat { 

/** 
 * A utility to check that the input `Rcpp::StringVector` is of length 1 and then to convert it into a `std::string`.
 *
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @param str An R object containing a character vector of length 1.
 * @return The first (and only) string in `str`, but as a `std::string`.
 * 
 * @internal
 */
inline std::string make_to_string(const Rcpp::RObject& str) {
    Rcpp::StringVector as_str(str);
    if (as_str.size()!=1) { 
        throw std::runtime_error("input RObject should contain a single string");
    }
    return Rcpp::as<std::string>(as_str[0]);
}

/**
 * Extract the class infomation as the `class` attribute of an S4 object.
 *
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @param incoming An R object, expected to be an instance of an S4 class.
 *
 * @return Another R object containing the value of the `class` attribute.
 * This should be interpretable as a character vector of length 1,
 * itself containing the `package` attribute to specify the package of origin.
 */
inline Rcpp::RObject get_class_object(const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        throw std::runtime_error("object has no 'class' attribute");
    }
    return incoming.attr("class");
}

/**
 * Extract the class name for an S4 object.
 *
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @param incoming An R object, expected to be an instance of an S4 class.
 *
 * @return The name of the class.
 * 
 * @internal
 */
inline std::string get_class_name(const Rcpp::RObject& incoming) {
    return make_to_string(get_class_object(incoming));
}

/**
 * Extract the package of origin for a given class.
 *
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @param incoming An R object containing S4 class information, typically the output of `get_class_object`.
 *
 * @return The name of the package of origin.
 */
inline std::string extract_class_package(const Rcpp::RObject& classname) {
    if (!classname.hasAttribute("package")) {
        throw std::runtime_error("class name has no 'package' attribute");
    }
    return make_to_string(classname.attr("package"));
}

/**
 * Extract the class name and its the package of origin for an instance of an S4 object.
 *
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @param incoming An R object, expected to be an instance of an S4 class.
 *
 * @return A `std::pair` containing the name of the class and the package of origin for the class.
 */
inline std::pair<std::string, std::string> get_class_package(const Rcpp::RObject& incoming) {
    Rcpp::RObject classname=get_class_object(incoming);
    return std::make_pair(make_to_string(classname), extract_class_package(classname));
}

/**
 * Translate SEXP numbers to plain-English types, for the supported integer, logical, character and double-precision types.
 * 
 * @note This is an internal function and should not be called directly by **beachmat** users.
 *
 * @param sexp_type The code number for a given SEXP type.
 *
 * @return The plain-English name for that type.
 */
inline std::string translate_type(int sexp_type) {
    std::string should_be;
    switch(sexp_type) {
        case REALSXP:
            should_be="double";
            break;
        case INTSXP:
            should_be="integer";
            break;
        case LGLSXP:
            should_be="logical";
            break;
        case STRSXP:
            should_be="character";
            break;
        default:
            std::stringstream err;
            err << "unsupported sexptype '" << sexp_type << "'";
            throw std::runtime_error(err.str());
    }
    return should_be;
}

}

#endif
