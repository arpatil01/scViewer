subroutine swap(j,k1,k2,shdswp,nadj,madj,x,y,ntot,eps,incAdj)

# The segment k1->k2 is a diagonal of a quadrilateral
# with a vertex at j (the point being added to the
# triangulation).  If the LOP is not satisfied, swap
# it for the other diagonal.
# Called by addpt.

implicit double precision(a-h,o-z)
dimension nadj(-3:ntot,0:madj), x(-3:ntot), y(-3:ntot)
logical shdswp


# If vertices k1 and k2 are not connected there is no diagonal to swap.
# This could happen if vertices j, k1, and k2 were colinear, but shouldn't.
call adjchk(k1,k2,shdswp,nadj,madj,ntot)
if(!shdswp) return

# Get the other vertex of the quadrilateral.
call pred(k,k1,k2,nadj,madj,ntot)  # If these aren't the same, then
call succ(kk,k2,k1,nadj,madj,ntot) # there is no other vertex.
if(kk!=k) {
    shdswp = .false.
    return
}

# Check whether the LOP is satisified; i.e. whether
# vertex k is outside the circumcircle of vertices j, k1, and k2
call qtest(j,k1,k,k2,shdswp,x,y,ntot,eps)

# Do the actual swapping.
if(shdswp) {
    call delet(k1,k2,nadj,madj,ntot)
    call insrt(j,k,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
}
return
end
