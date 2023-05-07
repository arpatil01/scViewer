#ifndef BEACHTEST_GET_ERRORS_H
#define BEACHTEST_GET_ERRORS_H
#include "beachtest.h"

/* This function tests the edge cases and error triggers. */

template <class T, class M>  
void get_errors(M ptr, const Rcpp::IntegerVector& mode) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];

    T stuff;
    if (Mode==0) {
        ptr->get_row(0, stuff.begin(), 0, 0); // Should not break.
        ptr->get_col(0, stuff.begin(), 0, 0); 
    } else if (Mode==1) {
        ptr->get_row(-1, stuff.begin()); // break!
    } else if (Mode==-1) {
        ptr->get_col(-1, stuff.begin()); // break!
    } else if (Mode==2) {
        ptr->get_row(0, stuff.begin(), 1, 0); // break!
    } else if (Mode==-2) {
        ptr->get_col(0, stuff.begin(), 1, 0); // break!
    } else if (Mode==3) {
        ptr->get_row(0, stuff.begin(), 0, -1); // break!
    } else if (Mode==-3) {
        ptr->get_col(0, stuff.begin(), 0, -1); // break!
    }
    return;
}

template <class T, class M>  
void get_multi_errors(M ptr, const Rcpp::IntegerVector& mode) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];

    T stuff;
    if (Mode==1 || Mode==-1) {
        Rcpp::IntegerVector thingy=Rcpp::IntegerVector::create(1, 0, 2);
        if (Mode > 0) {
            ptr->get_rows(thingy.begin(), thingy.size(), stuff.begin(), 0, 0); // break!
        } else {
            ptr->get_cols(thingy.begin(), thingy.size(), stuff.begin(), 0, 0); // break!
        }

    } else if (Mode==2 || Mode==-2) {
        Rcpp::IntegerVector thingy=Rcpp::IntegerVector::create(0, 1, -1);
        if (Mode > 0) { 
            ptr->get_rows(thingy.begin(), thingy.size(), stuff.begin(), 0, 0); // break!
        } else { 
            ptr->get_cols(thingy.begin(), thingy.size(), stuff.begin(), 0, 0); // break!
        }

    } else if (Mode==3 || Mode==-3) {
        Rcpp::IntegerVector thingy=Rcpp::IntegerVector::create(0, 1, 2);
        if (Mode > 0) { 
            ptr->get_rows(thingy.begin(), thingy.size(), stuff.begin(), 1, 0); // break!
        } else {
            ptr->get_cols(thingy.begin(), thingy.size(), stuff.begin(), 1, 0); // break!
        }

    } else if (Mode==4 || Mode==-4) {
        Rcpp::IntegerVector thingy=Rcpp::IntegerVector::create(0, 1, 2);
        if (Mode > 0) { 
            ptr->get_rows(thingy.begin(), thingy.size(), stuff.begin(), 0, -1); // break!
        } else {
            ptr->get_cols(thingy.begin(), thingy.size(), stuff.begin(), 0, -1); // break!
        }
    }
    return;
}

#endif
