## ----setup, echo=FALSE, message=FALSE, warning=FALSE--------------------------
options(digits = 4)
set.seed(1234)

## ----datgroundbeef, echo=TRUE-------------------------------------------------
library("fitdistrplus")
data("groundbeef")
str(groundbeef)

## ----figgroundbeef, fig.align='center', fig.width=7, fig.height=4, fig.cap="Histogram and CDF plots of an empirical distribution for a continuous variable (serving size from the `groundbeef` data set) as provided by the `plotdist` function."----
plotdist(groundbeef$serving, histo = TRUE, demp = TRUE)

## ----descgroundbeefplot, fig.align='center', fig.width=5, fig.height=5, fig.cap="Skewness-kurtosis plot for a continuous variable (serving size from the `groundbeef` data set) as provided by the `descdist` function."----
descdist(groundbeef$serving, boot = 1000)

## ----fitgroundbeef.weibull----------------------------------------------------
fw <- fitdist(groundbeef$serving, "weibull")
summary(fw)

## ----groundbeefcomp, fig.align='center', fig.width=7, fig.height=7, fig.cap="Four Goodness-of-fit plots for various distributions fitted to continuous data (Weibull, gamma and lognormal distributions fitted to serving sizes from the `groundbeef` data set) as provided by functions `denscomp`, `qqcomp`, `cdfcomp` and `ppcomp`."----
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
fg <- fitdist(groundbeef$serving, "gamma")
fln <- fitdist(groundbeef$serving, "lnorm")
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)

## ----fitendo, fig.align='center', fig.width=6, fig.height=6, fig.cap="CDF plot to compare the fit of four distributions to acute toxicity values of various organisms for the organochlorine pesticide endosulfan (`endosulfan` data set) as provided by the `cdfcomp` function, with CDF values in a logscale to emphasize discrepancies on the left tail."----
library(actuar)
data("endosulfan")
ATV <- endosulfan$ATV
fendo.ln <- fitdist(ATV, "lnorm")
fendo.ll <- fitdist(ATV, "llogis", start = list(shape = 1, scale = 500))
fendo.P <- fitdist(ATV, "pareto", start = list(shape = 1, scale = 500))
fendo.B <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1))
cdfcomp(list(fendo.ln, fendo.ll, fendo.P, fendo.B), xlogscale = TRUE, 
        ylogscale = TRUE, legendtext = c("lognormal", "loglogistic", "Pareto", "Burr"))

## ----quantilefitdist, echo=TRUE, fig=FALSE------------------------------------
quantile(fendo.B, probs = 0.05)
quantile(ATV, probs = 0.05)

## ----fendo.gof.print, echo=TRUE, fig=FALSE------------------------------------
gofstat(list(fendo.ln, fendo.ll, fendo.P, fendo.B), 
        fitnames = c("lnorm", "llogis", "Pareto", "Burr"))

## ----fitBurr.boot.echo, echo=TRUE---------------------------------------------
bendo.B <- bootdist(fendo.B, niter = 1001)
summary(bendo.B)

## ----bootstrap, fig.align='center', fig.width=6, fig.height=6, fig.cap="Bootstrappped values of parameters for a fit of the Burr distribution characterized by three parameters (example on the `endosulfan` data set) as provided by the plot of an object of class `bootdist`."----
plot(bendo.B)

## ----fitATV.lnorm.quantile, echo=TRUE-----------------------------------------
quantile(bendo.B, probs = 0.05)

## ----plotfitMGE, fig.align='center', fig.width=6, fig.height=6, fig.cap="Comparison of a lognormal distribution fitted by MLE and by MGE using two different goodness-of-fit distances: left-tail Anderson-Darling and left-tail Anderson Darling of second order (example with the `endosulfan` data set) as provided by the `cdfcomp` function, with CDF values in a logscale to emphasize discrepancies on the left tail."----
fendo.ln.ADL <- fitdist(ATV, "lnorm", method = "mge", gof = "ADL")
fendo.ln.AD2L <- fitdist(ATV, "lnorm", method = "mge", gof = "AD2L")
cdfcomp(list(fendo.ln, fendo.ln.ADL, fendo.ln.AD2L), 
  xlogscale = TRUE, ylogscale = TRUE, 
  main = "Fitting a lognormal distribution",
  xlegend = "bottomright", 
  legendtext = c("MLE", "Left-tail AD", "Left-tail AD 2nd order"))

## ----quantilefitdist2, echo=TRUE, fig=FALSE-----------------------------------
(HC5.estimates <- c(
  empirical = as.numeric(quantile(ATV, probs = 0.05)), 
  Burr = as.numeric(quantile(fendo.B, probs = 0.05)$quantiles), 
  lognormal_MLE = as.numeric(quantile(fendo.ln, probs = 0.05)$quantiles), 
  lognormal_AD2 = as.numeric(quantile(fendo.ln.ADL, probs = 0.05)$quantiles), 
  lognormal_AD2L = as.numeric(quantile(fendo.ln.AD2L, probs = 0.05)$quantiles)))

## ----danish.mme, echo=TRUE, eval=TRUE-----------------------------------------
data("danishuni")
str(danishuni)

## ----danishmme, fig.align='center', fig.width=7, fig.height=4, fig.cap="Comparison between MME and MLE when fitting a lognormal or a Pareto distribution to loss data from the `danishuni` data set."----
fdanish.ln.MLE <- fitdist(danishuni$Loss, "lnorm")
fdanish.ln.MME <- fitdist(danishuni$Loss, "lnorm", method = "mme", order = 1:2)
library(actuar)
fdanish.P.MLE <- fitdist(danishuni$Loss, "pareto", start = list(shape = 10, scale = 10), 
                         lower = 2+1e-6, upper = Inf)
memp <- function(x, order) sum(x^order) / length(x)
fdanish.P.MME <- fitdist(danishuni$Loss, "pareto", method = "mme", order = 1:2, memp = "memp", 
                         start = list(shape = 10, scale = 10), lower = c(2+1e-6, 2+1e-6), 
                         upper = c(Inf, Inf))

par(mfrow = c(1, 2))
cdfcomp(list(fdanish.ln.MLE, fdanish.ln.MME), legend = c("lognormal MLE", "lognormal MME"),
        main = "Fitting a lognormal distribution", xlogscale = TRUE, datapch = 20)
cdfcomp(list(fdanish.P.MLE, fdanish.P.MME), legend = c("Pareto MLE", "Pareto MME"), 
        main = "Fitting a Pareto distribution", xlogscale = TRUE, datapch = 20)

## ----danish.mme.pareto, echo=TRUE, fig=FALSE----------------------------------
gofstat(list(fdanish.ln.MLE, fdanish.P.MLE, fdanish.ln.MME, fdanish.P.MME), 
        fitnames = c("lnorm.mle", "Pareto.mle", "lnorm.mme", "Pareto.mme"))

## ----danishqme, fig.align='center', fig.width=6, fig.height=6, fig.cap="Comparison between QME and MLE when fitting a lognormal distribution to loss data from the `danishuni` data set."----
fdanish.ln.QME1 <- fitdist(danishuni$Loss, "lnorm", method = "qme", probs = c(1/3, 2/3))
fdanish.ln.QME2 <- fitdist(danishuni$Loss, "lnorm", method = "qme", probs = c(8/10, 9/10))
cdfcomp(list(fdanish.ln.MLE, fdanish.ln.QME1, fdanish.ln.QME2), 
        legend = c("MLE", "QME(1/3, 2/3)", "QME(8/10, 9/10)"), 
        main = "Fitting a lognormal distribution", xlogscale = TRUE, datapch = 20)

## ----optimmethod.gamma, echo=TRUE---------------------------------------------
data("groundbeef")
fNM <- fitdist(groundbeef$serving, "gamma", optim.method = "Nelder-Mead")
fBFGS <- fitdist(groundbeef$serving, "gamma", optim.method = "BFGS") 
fSANN <- fitdist(groundbeef$serving, "gamma", optim.method = "SANN")
fCG <- try(fitdist(groundbeef$serving, "gamma", optim.method = "CG",
                   control = list(maxit = 10000)))
if(inherits(fCG, "try-error")) {fCG <- list(estimate = NA)}

## ----optimmethod.customgenoud, echo=TRUE--------------------------------------
mygenoud <- function(fn, par, ...) 
{
   require(rgenoud)
   res <- genoud(fn, starting.values = par, ...)        
   standardres <- c(res, convergence = 0)
   return(standardres)
}

## ----optimmethod.customgenoud.fitdist, echo=TRUE, eval=TRUE-------------------
fgenoud <- mledist(groundbeef$serving, "gamma", custom.optim = mygenoud, nvars = 2, 
                   max.generations = 10, Domains = cbind(c(0, 0), c(10, 10)), 
                   boundary.enforcement = 1, hessian = TRUE, print.level = 0, P9 = 10)
cbind(NM = fNM$estimate, BFGS = fBFGS$estimate, SANN = fSANN$estimate, CG = fCG$estimate, 
      fgenoud = fgenoud$estimate)

## ----datsalinity, echo=TRUE---------------------------------------------------
data("salinity")
str(salinity)

## ----plotsalinity2, fig.align='center', fig.width=6, fig.height=6, fig.cap="Simple plot of censored raw data (72-hour acute salinity tolerance of riverine macro-invertebrates from the `salinity` data set) as ordered points and intervals."----
plotdistcens(salinity, NPMLE = FALSE)

## ----plotdistcens, echo=TRUE, fig=FALSE---------------------------------------
fsal.ln <- fitdistcens(salinity, "lnorm")
fsal.ll <- fitdistcens(salinity, "llogis", start = list(shape = 5, scale = 40))
summary(fsal.ln)
summary(fsal.ll)

## ----cdfcompcens, fig.align='center', fig.width=7, fig.height=7, fig.cap="Some goodness-of-fit plots for fits of a lognormal and a loglogistic distribution to censored data: LC50 values from the `salinity` data set."----
par(mfrow = c(2, 2))
cdfcompcens(list(fsal.ln, fsal.ll), legendtext = c("lognormal", "loglogistic "))
qqcompcens(fsal.ln, legendtext = "lognormal")
ppcompcens(fsal.ln, legendtext = "lognormal")
qqcompcens(list(fsal.ln, fsal.ll), legendtext = c("lognormal", "loglogistic "),
           main = "Q-Q plot with 2 dist.")

## ----dattoxocara, echo=TRUE---------------------------------------------------
data("toxocara")
str(toxocara)

## ----fittoxocara.poisnbinom, echo = TRUE, fig = FALSE-------------------------
(ftoxo.P <- fitdist(toxocara$number, "pois"))
(ftoxo.nb <- fitdist(toxocara$number, "nbinom"))

## ----fittoxocarapoisnbinom, fig.align='center', fig.width=7, fig.height=4, fig.cap="Comparison of the fits of a negative binomial and a Poisson distribution to numbers of *Toxocara cati* parasites from the `toxocara` data set."----
par(mfrow = c(1, 2))
denscomp(list(ftoxo.P, ftoxo.nb), legendtext = c("Poisson", "negative binomial"), fitlty = 1)
cdfcomp(list(ftoxo.P, ftoxo.nb), legendtext = c("Poisson", "negative binomial"), fitlty = 1)

## ----fittoxocara.poisnbinom.gof-----------------------------------------------
gofstat(list(ftoxo.P, ftoxo.nb), fitnames = c("Poisson", "negative binomial"))

