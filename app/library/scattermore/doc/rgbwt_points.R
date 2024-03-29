## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scattermore)
points <- cbind(rnorm(1e5), rnorm(1e5))

## -----------------------------------------------------------------------------
rgbwt <- scatter_points_rgbwt(points, RGBA= c(64,128,192,50), xlim=c(-5,5), ylim=c(-5,5))
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
plot(raster, interpolate=F)

## -----------------------------------------------------------------------------
blurred_rgbwt <- apply_kernel_rgbwt(rgbwt)
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(blurred_rgbwt))
plot(raster, interpolate=F)
gauss_blurred_rgbwt <- apply_kernel_rgbwt(rgbwt, filter="gauss")
gauss_raster <- rgba_int_to_raster(rgbwt_to_rgba_int(gauss_blurred_rgbwt))
plot(gauss_raster, interpolate=F)

## -----------------------------------------------------------------------------
v <- c(255, 0, 0, 100, 0, 255, 0, 25, 0, 0, 255, 50, 0, 0, 0, 100)
palette <- array(v, c(4, 4))
map <- rep(c(1,2,3,4), each=25000)
rgbwt <- scatter_points_rgbwt(points, map = map, palette = palette, xlim=c(-5,5), ylim=c(-5,5))
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
plot(raster, interpolate=F)

## -----------------------------------------------------------------------------
v = c(255, 0, 0, 100, 0, 255, 0, 10, 0, 0, 255, 10, 0, 0, 0, 0)
colors <- array(v, c(4, 100000))
rgbwt <- scatter_points_rgbwt(points, RGBA = colors)
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
plot(raster, interpolate=F)

