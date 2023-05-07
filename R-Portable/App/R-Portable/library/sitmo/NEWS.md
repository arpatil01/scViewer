# sitmo 2.0.2

## DEPLOYMENT

- Switched from using TravisCI to using GitHub Actions.

## BUGFIXES

- Addressed unbalanced chunk delimiters in vignette. (thanks [yihui/knitr#2057](https://github.com/yihui/knitr/issues/2057)!)
- Addressed compilation with Rcpp 1.0.7 introducing new stream handling. 

# sitmo 2.0.1

## CHANGES

- Modified `vandercorput.h` by adding a newline to the end of the file.
- Modified `Makevars{.win}` to use `$(SHLIB_OPENMP_CXXFLAGS)` instead of `$(SHLIB_OPENMP_CFLAGS)` in `PKG_CXXFLAGS`

## TESTING INFRASTRUCTURE

- Modified `.travis.yml` to compile using both cores instead of only one.

# sitmo 2.0.0

## NEW

- Added `threefry` and `vandercorput` (#9, #10).
- Added contributor Ralf Stubner (@rstub)

## CHANGES

- Modified `threefry` to be fully _C++11_ compliant. (#11, @rstub)
- Removed the `Rcpp.plugin.maker()`

# sitmo 1.2.2

## BUGFIXES

- Addressed import of `Rcpp.plugin.maker()` by using the _exported_ variable name,
  e.g. `::`, in place of the _internal_ variable name, e.g. `:::` (#7).
- Updated `sitmo` URL to point to `stdfin/random` (#8).

# sitmo 1.2.1

## BUGFIXES

- Removed extra ; in sitmo header to quiet compile warnings (#4, thanks @helske)
- Updated `sitmo_two_seeds()` src and documentation in "Deployment of `sitmo` within C++ Code"
  so that it uses the second seed for eng2 and returns an n x 2 matrix instead of n x 3. (#5, thanks @helske)

# sitmo 1.2.0

## CHANGES

- Added plugin registration for `Rcpp:::Rcpp.plugin.maker()` (#3)
- Added `CxxFlags()` and `sitmoCxxFlags()` functions to display `CXX_FLAGS`
  required by `sitmo`. (#3)
- Updated examples in README.Rmd and SITMO internal vignette to 
  use the Rcpp depends attribute. (#3)

## BUGFIXES

- Corrected a signed and unsigned integer comparison in 
  "Making a Uniform PRNG with `sitmo`" vignette.
- Fixed notation in "Making a Uniform PRNG with `sitmo`" vignette.


# sitmo 1.1.0

## CHANGES

- Added `src/init.c` to address R 3.4 C++ registration requirement (#2)
- Clarified content in "Making a Uniform PRNG with `sitmo`" vignette.

## BUG FIXES

- Addressed signed and unsigned integer comparison in sitmo header (#1)
- Corrected a URL that was problematic in "Deployment of `sitmo` within C++ Code" vignette.

# sitmo 1.0.0

## NEW

- Embedded `sitmo` header file in an R package.
- Provided code examples using `sitmo` header file.
- Released three vignettes detailing `sitmo` contents: 
    - Deployment of `sitmo` within C++ Code
    - Making a Uniform PRNG with `sitmo`
    - `sitmo`'s BigCrush Results
    