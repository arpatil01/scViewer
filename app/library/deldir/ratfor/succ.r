subroutine  succ(ksc,i,j,nadj,madj,ntot)

# Find the successor of j in the adjacency list of i.
# Called by addpt, initad, trifnd, swap, delout, dirseg, dirout.

implicit double precision(a-h,o-z)
dimension nadj(-3:ntot,0:madj)
dimension ndi(1)

# Set dummy integer for call to intpr(...).
ndi(1) = 0

n = nadj(i,0)

# If the adjacency list of i is empty, then clearly j has no successor
# in this adjacency list.  Something's wrong; stop.
if(n==0) {
    call intpr("Adjacency list of i is empty, and so cannot contain j.",-1,ndi,0)
    call rexit("Bailing out of succ.")
}

# The adjacency list of i is non-empty; search through it until j is found;
# add 1 to the location of j, and find the contents of this new location.
do k = 1,n {
    if(j==nadj(i,k)) {
        kp = k+1
        if(kp>n) kp = 1         # Take kp modulo n. (The adjacency list
        ksc = nadj(i,kp)        # is circular.)
        return
    }
}

# The adjacency list doesn't contain j.  Something's wrong.
ndi(1) = i
call intpr("i =",-1,ndi,1)
ndi(1) = j
call intpr("j =",-1,ndi,1)
call intpr("Adjacency list of i does not contain j.",-1,ndi,0)
call rexit("Bailing out of succ.")
end
