## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scattermore)

## -----------------------------------------------------------------------------
histogram <- scatter_points_histogram(cbind(rnorm(1e5), rnorm(1e5)), xlim=c(-5,5), ylim=c(-5,5))
image(histogram)

## -----------------------------------------------------------------------------
blurred_histogram <- apply_kernel_histogram(histogram, radius=4)
image(blurred_histogram)
gauss_blurred_histogram <- apply_kernel_histogram(histogram, filter="gauss", radius=4)
image(gauss_blurred_histogram)

## -----------------------------------------------------------------------------
rgbwt <- histogram_to_rgbwt(blurred_histogram)
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
plot(raster, interpolate=F)

## -----------------------------------------------------------------------------
v = c(255, 0, 0, 100, 0, 255, 0, 25, 0, 0, 255, 50, 0, 0, 0, 100)
palette = array(v, c(4, 4))

rgbwt <- histogram_to_rgbwt(blurred_histogram, RGBA=palette)
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
plot(raster, interpolate=F)

