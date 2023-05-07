# This script transplants the Annoy headers from RcppAnnoy into BiocNeighbors.
# Fundamentally, this is necessary because changes in the Annoy libaries (and
# thus RcppAnnoy versions) cause changes in the results. Thus, we want a
# constant version in BiocNeighbors during the lifetime of a single BioC
# release. See the discussion in jlmelville/uwot#69 for more details.

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)!=2) {
    stop("need <version> <target>")
}

version <- args[1]
if (version=="current") {
    options(repos = "http://cran.us.r-project.org")
    df <- as.data.frame(available.packages())
    df <- df["RcppAnnoy",]
    url <- sprintf("https://cran.r-project.org/src/contrib/RcppAnnoy_%s.tar.gz", df$Version)
} else {
    url <- sprintf("https://cran.r-project.org/src/contrib/Archive/RcppAnnoy/RcppAnnoy_%s.tar.gz", version)
}

tmp <- tempfile(fileext=".tar.gz")
download.file(url, tmp)
exdir <- tempfile()
untar(tmp, exdir=exdir)

target <- args[2]
dest <- file.path(target, "src/annoy")
unlink(dest, recursive=TRUE)
dir.create(dest)

src <- file.path(exdir, "RcppAnnoy/inst/include")
contents <- list.files(src)
for (x in contents) {
    file.copy(file.path(src, x), file.path(dest, x), recursive=TRUE)
}
