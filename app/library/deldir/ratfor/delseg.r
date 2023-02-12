subroutine delseg(delsgs,ndel,nadj,madj,nn,x,y,ntot,incSeg)

# Output the endpoints of the line segments joining the
# vertices of the Delaunay triangles.
# Called by master.

implicit double precision(a-h,o-z)
logical value
dimension nadj(-3:ntot,0:madj), x(-3:ntot), y(-3:ntot)
dimension delsgs(6,ndel)

# Initialise incSeg
incSeg = 0

# For each distinct pair of points i and j, if they are adjacent
# then put their endpoints into the output array.
nn = ntot-4
kseg = 0
do i = 2,nn {
    do j = 1,i-1 {
        call adjchk(i,j,value,nadj,madj,ntot)
        if(value) {
            kseg = kseg+1
            if(kseg > ndel) {
                incSeg = 1
                return
            }
            delsgs(1,kseg) = x(i)
            delsgs(2,kseg) = y(i)
            delsgs(3,kseg) = x(j)
            delsgs(4,kseg) = y(j)
            delsgs(5,kseg) = i
            delsgs(6,kseg) = j
        }
    }
}
ndel = kseg

return
end
