## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scattermore)
points <- cbind(rnorm(1e5), rnorm(1e5))
p1 <- scatter_points_rgbwt(points, RGBA= c(64,128,192,50), xlim=c(-5,5), ylim=c(-5,5))
p2 <- scatter_points_rgbwt(points, RGBA= c(192,128,64,50), xlim=c(-5,5), ylim=c(-5,5))

## -----------------------------------------------------------------------------
l <- list(p1, p2)
merged <- merge_rgbwt(l)
raster <- rgba_int_to_raster(rgbwt_to_rgba_int(merged))
plot(raster, interpolate=F)

## -----------------------------------------------------------------------------
p1_frgba <- rgbwt_to_rgba_float(p1)
p2_frgba <- rgbwt_to_rgba_float(p2)
l <- list(p1_frgba, p2_frgba)
blended <- blend_rgba_float(l)
raster <- rgba_int_to_raster(rgba_float_to_rgba_int(blended))
plot(raster, interpolate=F)

