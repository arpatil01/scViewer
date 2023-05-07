#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocGenerics rbind cbind t
delayed_funs <- function(basefun, DELAYED_FUN)
# Creates a list of supported DelayedMatrix operations 
# that are handled natively in beachmat, without the 
# need for explicit realization.
{
    # Subsetting:
    .choose_half <- function(n) {
        sample(n, ceiling(n/2))
    }

    subset_row_fun <- function(...) {
        out <- DelayedArray(basefun(...))
        out[.choose_half(nrow(out)),]
    }

    subset_col_fun <- function(...) {
        out <- DelayedArray(basefun(...))
        out[,.choose_half(ncol(out))]
    }

    subset_both_fun <- function(...) {
        out <- DelayedArray(basefun(...))
        out[.choose_half(nrow(out)),.choose_half(ncol(out))]
    }

    # Dimnaming
    .namers <- function(n, pre) {
        sprintf("%s%i", pre, seq_len(n))
    }

    name_row_fun <- function(...) { 
        out <- DelayedArray(basefun(...))
        rownames(out) <- .namers(nrow(out), "Gene")
        return(out)
    }

    name_col_fun <- function(...) { 
        out <- DelayedArray(basefun(...))
        colnames(out) <- .namers(ncol(out), "Cell")
        return(out)
    }

    name_both_fun <- function(...) {
        out <- DelayedArray(basefun(...))
        rownames(out) <- .namers(nrow(out), "Gene")
        colnames(out) <- .namers(ncol(out), "Cell")
        return(out)
    }

    # Transposition
    trans_fun <- function(...) { 
        out <- DelayedArray(basefun(...))
        DelayedArray::t(out)
    }

    trans_subset_fun <- function(...) {
        t(subset_both_fun(...))
    }

    # Combining to generate an odd seed type that switches to chunked access.
    rcomb_fun <- function(...) {
        rbind(DelayedArray(basefun(...)), DelayedArray(basefun(...)))
    }

    ccomb_fun <- function(...) {
        cbind(DelayedArray(basefun(...)), DelayedArray(basefun(...)))
    }

    # Performing a delayed operation (DELAYED_FUN needs to be chosen to ensure the modification is of the same type).
    dops_fun <- function(...) { 
        out <- DelayedArray(basefun(...))
        DELAYED_FUN(out)
    }

    dops_rsubset_fun <- function(...) {
        out <- subset_row_fun(...)
        DELAYED_FUN(out)
    }

    dops_csubset_fun <- function(...) {
        out <- subset_col_fun(...)
        DELAYED_FUN(out)
    }

    dops_transpose_fun <- function(...) { # Checking that transposition WITH delayed ops wipes the transformer.
        t(dops_fun(...))
    }

    return(list(SR=subset_row_fun,
                SC=subset_col_fun,
                SB=subset_both_fun,
                NR=name_row_fun,
                NC=name_col_fun,
                NB=name_both_fun,
                TR=trans_fun,
                TRS=trans_subset_fun,
                RC=rcomb_fun,
                CC=ccomb_fun,
                DO=dops_fun,
                DOR=dops_rsubset_fun,
                DOC=dops_csubset_fun,
                DOT=dops_transpose_fun))
}

