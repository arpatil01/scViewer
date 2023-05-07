#ifndef BEACHTEST_SET_ERRORS_H
#define BEACHTEST_SET_ERRORS_H
#include "beachtest.h"

/* This function tests the edge cases and error triggers. */

template <class T, class M>  
void set_errors(M ptr, const Rcpp::IntegerVector& mode) {
    if (mode.size()!=1) { 
        throw std::runtime_error("'mode' should be an integer scalar"); 
    }
    const int Mode=mode[0];

    T stuff;
    if (Mode==0) {
        ptr->set_row(0, stuff.begin(), 0, 0); // Should not break.
        ptr->set_col(0, stuff.begin(), 0, 0); 
    } else if (Mode==1) {
        ptr->set_row(-1, stuff.begin(), 0, 0); // break!
    } else if (Mode==-1) {
        ptr->set_col(-1, stuff.begin(), 0, 0); // break!
    } else if (Mode==2) {
        ptr->set_row(0, stuff.begin(), 1, 0); // break!
    } else if (Mode==-2) {
        ptr->set_col(0, stuff.begin(), 1, 0); // break!
    } else if (Mode==3) {
        ptr->set_row(0, stuff.begin(), 0, -1); // break!
    } else if (Mode==-3) {
        ptr->set_col(0, stuff.begin(), 0, -1); // break!
    }
   
    return;
}

#endif
