subroutine delet(i,j,nadj,madj,ntot)
# Delete i and j from each other's adjacency lists.
# Called by initad, swap.

implicit double precision(a-h,o-z)
dimension nadj(-3:ntot,0:madj)
logical adj

# First check that they're IN each other's lists.
call adjchk(i,j,adj,nadj,madj,ntot)

# Then do the actual deletion if they are.
if(adj) {
    call delet1(i,j,nadj,madj,ntot)
    call delet1(j,i,nadj,madj,ntot)
}

return
end
