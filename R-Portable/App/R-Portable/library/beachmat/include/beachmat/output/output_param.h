#ifndef BEACHMAT_OUTPUT_PARAM_H
#define BEACHMAT_OUTPUT_PARAM_H

#include "Rcpp.h"

#include "../utils/utils.h"

namespace beachmat {

class output_param {
public:
    output_param() {} // TODO: add option to fetch getRealizationBackend().

    output_param(const std::string& m, const std::string& p) : cls(m), pkg(p) {}

    output_param(Rcpp::RObject in) {
        if (!in.isS4()) {
            return;
        }
        auto classinfo=get_class_package(in);
        cls=classinfo.first;
        pkg=classinfo.second;
        return;
    }

    output_param(SEXP in) : output_param(Rcpp::RObject(in)) {}

    template<class M>
    explicit output_param(M mat) : output_param(mat->get_class(), mat->get_package()) {}

    bool is_external_available(const std::string& type) const {
        // Skipping packages that we know won't support this.
        if (pkg=="" || pkg=="Matrix" || pkg=="base" || pkg=="DelayedArray") { return false; }
        return has_external_support(type, cls, pkg, "output");
    }

    std::string get_class() const { return cls; }
    std::string get_package() const { return pkg; }
private:
    std::string cls="matrix";
    std::string pkg="base";
};

}

#endif
