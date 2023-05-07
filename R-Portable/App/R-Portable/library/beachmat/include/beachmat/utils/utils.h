#ifndef BEACHMAT_UTILS_H
#define BEACHMAT_UTILS_H

#include "Rcpp.h"

#include <string>
#include <sstream>
#include <utility>
#include <tuple>
#include <stdexcept>

namespace beachmat { 

template<class V>
using const_col_indexed_info=std::tuple<size_t, Rcpp::IntegerVector::iterator, typename V::iterator>;

/* String-related helper functions */

inline std::string make_to_string(const Rcpp::RObject& str) {
    Rcpp::StringVector as_str(str);
    if (as_str.size()!=1) { 
        throw std::runtime_error("input RObject should contain a single string");
    }
    return Rcpp::as<std::string>(as_str[0]);
}

/* Class checks. */

inline Rcpp::RObject get_class_object(const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        throw std::runtime_error("object has no 'class' attribute");
    }
    return incoming.attr("class");
}

inline std::string get_class_name(const Rcpp::RObject& incoming) {
    return make_to_string(get_class_object(incoming));
}

inline std::string extract_class_package(const Rcpp::RObject& classname) {
    if (!classname.hasAttribute("package")) {
        throw std::runtime_error("class name has no 'package' attribute");
    }
    return make_to_string(classname.attr("package"));
}

inline std::pair<std::string, std::string> get_class_package(const Rcpp::RObject& incoming) {
    Rcpp::RObject classname=get_class_object(incoming);
    return std::make_pair(make_to_string(classname), extract_class_package(classname));
}

inline Rcpp::RObject get_safe_slot(const Rcpp::RObject& incoming, const std::string& slotname) {
    if (!incoming.hasSlot(slotname)) { 
        throw std::runtime_error(std::string("no '") + slotname + "' slot in the " + get_class_name(incoming) + " object");
    }
    return incoming.slot(slotname);
}

inline void quit_on_df (const std::string& classname) {
    if (classname=="data.frame") {
        throw std::runtime_error("data.frames should be converted to matrices");
    }
    return;
}

inline void quit_on_df (const Rcpp::RObject& incoming) {
    if (incoming.isObject()) {
        quit_on_df(get_class_name(incoming));
    }
    return;
}

/* Type checks */

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

inline int find_sexp_type (const Rcpp::RObject& incoming) {
    if (!incoming.isObject()) {
        return incoming.sexp_type();
    }
    
    const auto classinfo=get_class_object(incoming);
    const std::string classname=make_to_string(classinfo);
    quit_on_df(classname);

    if (extract_class_package(classinfo)=="Matrix" && classname.length()==9 && classname.substr(3)=="Matrix") {
        if (classname[0]=='d') {
            return REALSXP;
        } else if (classname[0]=='l') {
            return LGLSXP;
        }

    } else {
        Rcpp::Environment delayenv=Rcpp::Environment::namespace_env("BiocGenerics");
        Rcpp::Function typefun=delayenv["type"];
        std::string curtype=Rcpp::as<std::string>(typefun(incoming));
        if (curtype=="logical") {
            return LGLSXP;
        } else if (curtype=="character") {
            return STRSXP;
        } else if (curtype=="integer") {
            return INTSXP;
        } else if (curtype=="double") {
            return REALSXP;
        }
    } 
    throw std::runtime_error(std::string("unknown SEXP type for ") + classname + " object");
}

}

#endif
