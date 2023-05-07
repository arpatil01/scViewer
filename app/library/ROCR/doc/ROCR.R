## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ROCR)

## ---- echo = FALSE, results = 'asis'------------------------------------------
table <- data.frame(group = c("Contingency ratios",
                              "Discrete covariation measures",
                              "Information retrieval measures",
                              "Performance in ROC space",
                              "Absolute scoring performance",
                              "Cost measures"),
                    measure = c("error rate, accuracy, sensitivity, specificity, true/false positive rate, fallout, miss, precision, recall, negative predictive value, prediction-conditioned fallout/miss.",
                                "Phi/Matthews correlation coefficient, mutual information, Chi-squared test statistic, odds ratio",
                                "F-measure, lift, precision-recall break-even point",
                                "ROC convex hull, area under the ROC curve",
                                "calibration error, mean cross-entropy, root mean-squared error",
                                "expected cost, explicit cost"))
knitr::kable(table,
             caption = "***Table 1:**Performance measures in the ROCR package*",
             col.names = c("",""),
             align = "l")

## -----------------------------------------------------------------------------
data(ROCR.hiv)
predictions <- ROCR.hiv$hiv.svm$predictions
labels <- ROCR.hiv$hiv.svm$labels
pred <- prediction(predictions, labels)
pred

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
perf <- performance(pred, "tpr", "fpr")
perf
plot(perf,
     avg="threshold",
     spread.estimate="boxplot")

## ---- echo=FALSE, results='asis', fig.asp=0.35, fig.width=7, fig.align='center',fig.cap="***Fig 1:** Visualizations of classifier performance (HIV coreceptor usage data): (a) receiver operating characteristic (ROC) curve; (b) peak accuracy across a range of cutoffs; (c) absolute difference between empirical and predicted rate of positives for windowed cutoff ranges, in order to evaluate how well the scores are calibrated as probability estimates. Owing to the probabilistic interpretation, cutoffs need to be in the interval [0,1], in contrast to other performance plots. (d) Score density estimates for the negative (solid) and positive (dotted) class.*"----
data(ROCR.hiv)
pp.unnorm <- ROCR.hiv$hiv.svm$predictions
ll <- ROCR.hiv$hiv.svm$labels

# normalize predictions to 0..1
v <- unlist(pp.unnorm)
pp <- lapply(pp.unnorm, function(run) {approxfun(c(min(v), max(v)), c(0,1))(run)})

par(mfrow=c(1,4))
pred<- prediction(pp, ll)
perf <- performance(pred, "tpr", "fpr")

plot(perf, avg= "threshold", colorize=TRUE, lwd= 3,
     coloraxis.at=seq(0,1,by=0.2),)
plot(perf, col="gray78", add=TRUE)
plot(perf, avg= "threshold", colorize=TRUE, colorkey=FALSE,lwd= 3,,add=TRUE)
mtext(paste0("(a)"), side = 3, adj = 0.01,line = 1)

perf <- performance(pred, "acc")
plot(perf, avg= "vertical", spread.estimate="boxplot", lwd=3,col='blue',
     show.spread.at= seq(0.1, 0.9, by=0.1),)
mtext(paste0("(b)"), side = 3, adj = 0.01,line = 1)


plot(performance(pred, "cal", window.size= 10),
     avg="vertical",)
mtext(paste0("(c)"), side = 3, adj = 0.01,line = 1)

plot(0,0,type="n", xlim= c(0,1), ylim=c(0,7),
     xlab="Cutoff", ylab="Density",)
mtext(paste0("(d)"), side = 3, adj = 0.01,line = 1)
for (runi in 1:length(pred@predictions)) {
  lines(density(pred@predictions[[runi]][pred@labels[[runi]]=="-1"]), col= "red")
  lines(density(pred@predictions[[runi]][pred@labels[[runi]]=="1"]), col="green")
}

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
perf <- performance(pred, "tpr", "fpr")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "With ROCR you can produce standard plots\nlike ROC curves ...")
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
perf <- performance(pred, "prec", "rec")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "... Precision/Recall graphs ...")
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
perf <- performance(pred, "sens", "spec")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main="... Sensitivity/Specificity plots ...")
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
perf <- performance(pred, "lift", "rpp")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "... and Lift charts.")
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)

## -----------------------------------------------------------------------------
data(ROCR.xval)
predictions <- ROCR.xval$predictions
labels <- ROCR.xval$labels
length(predictions)

## -----------------------------------------------------------------------------
pred <- prediction(predictions, labels)
perf <- performance(pred,'tpr','fpr')

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
plot(perf,
     colorize=TRUE,
     lwd=2,
     main='ROC curves from 10-fold cross-validation')

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
plot(perf,
     avg='vertical',
     spread.estimate='stderror',
     lwd=3,main='Vertical averaging + 1 standard error',
     col='blue')

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
plot(perf,
     avg='horizontal',
     spread.estimate='boxplot',
     lwd=3,
     main='Horizontal averaging + boxplots',
     col='blue')

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
plot(perf,
     avg='threshold',
     spread.estimate='stddev',
     lwd=2,
     main='Threshold averaging + 1 standard deviation',
     colorize=TRUE)

## ---- fig.asp=1, fig.width=6, fig.align='center'------------------------------
plot(perf,
     print.cutoffs.at=seq(0,1,by=0.2),
     text.cex=0.8,
     text.y=lapply(as.list(seq(0,0.5,by=0.05)), function(x) { rep(x,length(perf@x.values[[1]])) } ),
     col= as.list(terrain.colors(10)),
     text.col= as.list(terrain.colors(10)), 
     points.col= as.list(terrain.colors(10)), 
     main= "Cutoff stability")

## -----------------------------------------------------------------------------
perf <- performance(pred,"pcmiss","lift")

## ---- fig.asp=1, fig.width=5, fig.align='center'------------------------------
plot(perf,
     colorize=TRUE,
     print.cutoffs.at=seq(0,1,by=0.1),
     text.adj=c(1.2,1.2),
     avg="threshold",
     lwd=3,
     main= "You can freely combine performance measures ...")

