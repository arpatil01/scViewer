subroutine binsrt(x,y,rw,nn,ind,rind,tx,ty,ilst)
# Sort the data points into bins.
# Called by master.

implicit double precision(a-h,o-z)
dimension x(nn), y(nn), tx(nn), ty(nn)
integer rind(nn)
dimension ind(nn), ilst(nn)
dimension rw(4)
dimension ndi(1)

# Set dummy integer for call to intpr(...).
ndi(1) = 0

kdiv   = int(1+dble(nn)**0.25) # Round high.
xkdiv  = dble(kdiv)

# Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

w = xmax-xmin
h = ymax-ymin

# Number of bins is to be approx. sqrt(nn); thus number of subdivisions
# on each side of rectangle is approx. nn**(1/4).
dw  = w/xkdiv
dh  = h/xkdiv

# The width of each bin is dw; the height is dh.  We shall move across
# the rectangle from left to right, then up, then back from right to
# left, then up, ....  Note that kx counts the divisions from the left,
# ky counts the divisions from the bottom; kx is incremented by ink, which
# is +/- 1 and switches sign when ky is incremented; ky is always
# incremented by 1.
kx   = 1
ky   = 1
ink  = 1
k    = 0
do i = 1,nn { ilst(i) = 0 } # Keeps a list of those points already added
while(ky<=kdiv) {            # to the new list.
    do i = 1,nn {
        if(ilst(i)==1) next  # The i-th point has already been added
                             # to the new list.
# If the i-th point is in the current bin, add it to the list.
        xt = x(i)
        yt = y(i)
        ix = int(1+(xt-xmin)/dw)
        if(ix>kdiv) ix = kdiv
        jy = int(1+(yt-ymin)/dh)
        if(jy>kdiv) jy = kdiv
        if(ix==kx&jy==ky) {
            k = k+1
            ind(i)  = k  # Index i is the pos'n. of (x,y) in the
            rind(k) = i  # old list; k is its pos'n. in the new one.
            tx(k)   = xt
            ty(k)   = yt
            ilst(i) = 1  # Cross the i-th point off the old list.
        }
    }
# Move to the next bin.
    kc = kx+ink
    if((1<=kc)&(kc<=kdiv)) kx = kc
    else {
        ky  = ky+1
        ink = -ink
    }
}

# Check that all points from old list have been added to the new,
# with no spurious additions.
if(k!=nn) {
    call intpr("Mismatch between number of points",-1,ndi,0)
    call intpr("and number of sorted points.",-1,ndi,0)
    call rexit("Bailing out of binsrt.")
}

# Copy the new sorted vectors back on top of the old ones.
do i = 1,nn {
    x(i)   = tx(i)
    y(i)   = ty(i)
}

return
end
