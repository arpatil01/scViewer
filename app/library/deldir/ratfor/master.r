subroutine master(x,y,rw,nn,ntot,nadj,madj,eps,delsgs,ndel,delsum,
                  dirsgs,ndir,dirsum,incAdj,incSeg)

# Master subroutine:
# One subroutine to rule them all,
# One subroutine to find them.
# One subroutine to bring them all in,
# And in the darkness bind them.

# Note: "incAdj" <--> increase size of adjacency list.
#       "incSeg" <--> increase size of storage for segments.

implicit double precision(a-h,o-z)
dimension x(-3:ntot), y(-3:ntot)
dimension nadj(-3:ntot,0:madj)
dimension rw(4)
dimension delsgs(6,ndel), dirsgs(10,ndir)
dimension delsum(nn,4), dirsum(nn,3)

# Define one.
one = 1.d0

# Initialize the adjacency list; counts to 0, other entries to -99.
do i = -3,ntot {
    nadj(i,0) = 0
    do j = 1,madj {
        nadj(i,j) = -99
    }
}

# Put the four ideal points into x and y and the adjacency list.
# The ideal points are given pseudo-coordinates
# (-1,-1), (1,-1), (1,1), and (-1,1).  They are numbered as
#    0       -1      -2         -3
# i.e. the numbers decrease anticlockwise from the
# `bottom left corner'.
x(-3) = -one
y(-3) =  one
x(-2) =  one
y(-2) =  one
x(-1) =  one
y(-1) = -one
x(0)  = -one
y(0)  = -one

do i = 1,4 {
    j = i-4
    k = j+1
    if(k>0) k = -3
    call insrt(j,k,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
}

# Put in the first of the point set into the adjacency list.
do i = 1,4 {
    j = i-4
    call insrt(1,j,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
}
ntri = 4

# Now add the rest of the point set
do j = 2,nn {
    call addpt(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)
    if(incAdj==1) return
    ntri = ntri + 3
}

# Obtain the description of the triangulation.
call delseg(delsgs,ndel,nadj,madj,nn,x,y,ntot,incSeg)
if(incSeg==1) return

call delout(delsum,nadj,madj,x,y,ntot,nn)

call dirseg(dirsgs,ndir,nadj,madj,nn,x,y,ntot,rw,eps,ntri,incAdj,incSeg)
if(incAdj==1 | incSeg==1) return
call dirout(dirsum,nadj,madj,x,y,ntot,nn,rw,eps)
return
end
