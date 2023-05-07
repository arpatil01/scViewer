## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scattermore)
lines <- matrix(rnorm(40000),ncol=4,byrow=F)

## -----------------------------------------------------------------------------
rgbwt <- scatter_lines_rgbwt(lines, RGBA= c(64,128,192,50), xlim=c(-5,5), ylim=c(-5,5))
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
plot(raster, interpolate=F)

## -----------------------------------------------------------------------------
blurred_rgbwt <- apply_kernel_rgbwt(rgbwt)
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(blurred_rgbwt))
plot(raster, interpolate=F)
gauss_blurred_rgbwt <- apply_kernel_rgbwt(rgbwt, filter="gauss")
gauss_raster <- rgba_int_to_raster(rgbwt_to_rgba_int(gauss_blurred_rgbwt))
plot(gauss_raster, interpolate=F)

