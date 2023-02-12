subroutine dirseg(dirsgs,ndir,nadj,madj,nn,x,y,ntot,rw,eps,ntri,incAdj,incSeg)

# Output the endpoints of the segments of boundaries of Dirichlet
# tiles.  (Do it economically; each such segment once and only once.)
# Called by master.

implicit double precision(a-h,o-z)
logical collin, adjace, intfnd, bptab, bptcd, goferit, rwu
dimension nadj(-3:ntot,0:madj), x(-3:ntot), y(-3:ntot)
dimension dirsgs(10,ndir), rw(4)
dimension ndi(1)

# Set dummy integer for call to intpr(...).
ndi(1) = 0

# Initialise incSeg
incSeg = 0

# Add in some dummy corner points, outside the actual window.
# Far enough out so that no resulting tile boundaries intersect the
# window.

# Note that these dummy corners are needed by the routine `dirout'
# but will screw things up for `delseg' and `delout'.  Therefore
# this routine (`dirseg') must be called ***before*** dirout, and
# ***after*** delseg and delout.

# Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

a = xmax-xmin
b = ymax-ymin
c = sqrt(a*a+b*b)

nn  = ntot-4
nstt = nn+1
i = nstt
x(i) = xmin-c
y(i) = ymin-c
i = i+1
x(i) = xmax+c
y(i) = ymin-c
i = i+1
x(i) = xmax+c
y(i) = ymax+c
i = i+1
x(i) = xmin-c
y(i) = ymax+c

do j = nstt,ntot {
    call addpt(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)
    if(incAdj==1) return
    ntri = ntri + 3
}

# Put the segments into the array dirsgs.

# For each distinct pair of (genuine) data points, find out if they are
# adjacent.  If so, find the circumcentres of the triangles lying on each
# side of the segment joining them.
kseg = 0
do i = 2,nn {
    do j = 1,i-1 {
        call adjchk(i,j,adjace,nadj,madj,ntot)
        if(adjace) {
            call pred(k,i,j,nadj,madj,ntot)
            call circen(i,k,j,a,b,x,y,ntot,eps,collin)
            if(collin) {
                call intpr("Vertices of triangle are collinear.",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            }
            call succ(l,i,j,nadj,madj,ntot)
            call circen(i,j,l,c,d,x,y,ntot,eps,collin)
            if(collin) {
                call intpr("Vertices of triangle are collinear.",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            }
# If a circumcentre is outside the rectangular window
# of interest, draw a line joining it to the other
# circumcentre.  Find the intersection of this line with
# the boundary of the window; for (a,b) and call the point
# of intersection (ai,bi).  For (c,d), call it (ci,di).
# Note: rwu = "right way up".
            xi = x(i)
            xj = x(j)
            yi = y(i)
            yj = y(j)
            if(yi!=yj) {
                slope = (xi - xj)/(yj - yi)
                rwu   = .true.
            } else {
                slope = 0.d0
                rwu = .false.
            }
            call dldins(a,b,slope,rwu,ai,bi,rw,intfnd,bptab,nedgeab)
            if(!intfnd) {
                call intpr("Line from midpoint to circumcenter",-1,ndi,0)
                call intpr("does not intersect rectangle boundary!",-1,ndi,0)
                call intpr("But it HAS to!!!",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            }
            call dldins(c,d,slope,rwu,ci,di,rw,intfnd,bptcd,nedgecd)
            if(!intfnd) {
                call intpr("Line from midpoint to circumcenter",-1,ndi,0)
                call intpr("does not intersect rectangle boundary!",-1,ndi,0)
                call intpr("But it HAS to!!!",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            }
            goferit = .false.
            if(bptab & bptcd) {
                xm = 0.5*(ai+ci)
                ym = 0.5*(bi+di)
                if(xmin<xm&xm<xmax&ymin<ym&ym<ymax) {
                    goferit = .true.
                }
            }
            if((!bptab)|(!bptcd)) goferit = .true.
            if(goferit) {
                kseg = kseg + 1
                if(kseg > ndir) {
                    incSeg = 1
                    return
                }
                dirsgs(1,kseg) = ai
                dirsgs(2,kseg) = bi
                dirsgs(3,kseg) = ci
                dirsgs(4,kseg) = di
                dirsgs(5,kseg) = i
                dirsgs(6,kseg) = j
                if(bptab) dirsgs(7,kseg) = 1.d0
                else dirsgs(7,kseg) = 0.d0
                if(bptcd) dirsgs(8,kseg) = 1.d0
                else dirsgs(8,kseg) = 0.d0
                if(bptab) dirsgs(9,kseg) = -nedgeab
                else dirsgs(9,kseg) = k
                if(bptcd) dirsgs(10,kseg) = -nedgecd
                else dirsgs(10,kseg) = l
            }
        }
    }
}
ndir = kseg

return
end
