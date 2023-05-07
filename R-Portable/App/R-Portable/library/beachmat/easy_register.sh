pkgname=$1
shift
tmp="$*"
matrices=$(echo $tmp | sed "s/ /|/g")
src=$(ls *.cpp | grep -v "exports.cpp")

# Creating the export header.
cat << EOT > exports.h
#ifndef EXPORTS_H
#define EXPORTS_H
#include "Rcpp.h"

extern "C" {

EOT

cat ${src} | egrep "(${matrices})_[^_]+_(input|output)_.*{" | sed -E "s/ [^ ]+,/,/g" | sed -E "s/ [^ ]+\) *\{/\);/" | sed "s/;/;\n/" >> exports.h

cat << EOT >> exports.h
}

#endif
EOT

# Creating the export file.
cat << EOT > exports.cpp
#include "exports.h"
#include "R_ext/Rdynload.h"

#define REGISTER(x) R_RegisterCCallable("$pkgname", #x, reinterpret_cast<DL_FUNC>(x))

extern "C" {

void R_init_$pkgname(DllInfo *info) {

EOT

cat ${src} | egrep "(${matrices})_[^_]+_(input|output)_.*{" | sed -E "s/ ?\(.*$//g" | sed -E "s/^.* //" | sed "s/^\(.*\)$/REGISTER(\1);\n/" >> exports.cpp

cat << EOT >> exports.cpp
}

}
EOT
