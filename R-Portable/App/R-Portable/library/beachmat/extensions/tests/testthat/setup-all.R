testpkg <- system.file("testpkg", package="beachmat")
devtools::install(testpkg, build=FALSE, upgrade="never", reload=FALSE, dependencies=FALSE)
library(beachtest)
