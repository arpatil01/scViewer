## ---- echo=FALSE, results="hide"----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE)
library(BiocStyle)
set.seed(10918)

## ----quickstart-load-data, message=FALSE, warning=FALSE-----------------------
library(scRNAseq)
example_sce <- ZeiselBrainData()
example_sce

## -----------------------------------------------------------------------------
library(scater)
example_sce <- addPerCellQC(example_sce, 
    subsets=list(Mito=grep("mt-", rownames(example_sce))))

## -----------------------------------------------------------------------------
plotColData(example_sce, x = "sum", y="detected", colour_by="tissue") 

## ----plot-pdata-pct-exprs-controls--------------------------------------------
plotColData(example_sce, x = "sum", y="subsets_Mito_percent", 
    other_fields="tissue") + facet_wrap(~tissue)

## ----plot-highest, fig.asp=1, fig.wide=TRUE, eval=.Platform$OS.type!="windows"----
#  plotHighestExprs(example_sce, exprs_values = "counts")

## -----------------------------------------------------------------------------
# Computing variance explained on the log-counts,
# so that the statistics reflect changes in relative expression.
example_sce <- logNormCounts(example_sce)

vars <- getVarianceExplained(example_sce,
    variables=c("tissue", "total mRNA mol", "sex", "age"))
head(vars)

plotExplanatoryVariables(vars)

## ----plot-expression, fig.wide=TRUE-------------------------------------------
plotExpression(example_sce, rownames(example_sce)[1:6], x = "level1class")

## ----plot-expression-scatter--------------------------------------------------
plotExpression(example_sce, rownames(example_sce)[1:6],
    x = rownames(example_sce)[10])

## ----plot-expression-col------------------------------------------------------
plotExpression(example_sce, rownames(example_sce)[1:6],
    x = "level1class", colour_by="tissue")

## ----plot-expression-many-----------------------------------------------------
plotExpression(example_sce, rownames(example_sce)[1:6])

## -----------------------------------------------------------------------------
example_sce <- runPCA(example_sce)
str(reducedDim(example_sce, "PCA"))

## -----------------------------------------------------------------------------
example_sce <- runPCA(example_sce, name="PCA2",
    subset_row=rownames(example_sce)[1:1000],
    ncomponents=25)
str(reducedDim(example_sce, "PCA2"))

## ----plot-tsne-1comp-colby-sizeby-exprs---------------------------------------
# Perplexity of 10 just chosen here arbitrarily.
set.seed(1000)
example_sce <- runTSNE(example_sce, perplexity=10)
head(reducedDim(example_sce, "TSNE"))

## ----plot-tsne-from-pca-------------------------------------------------------
set.seed(1000)
example_sce <- runTSNE(example_sce, perplexity=50, 
    dimred="PCA", n_dimred=10)
head(reducedDim(example_sce, "TSNE"))

## -----------------------------------------------------------------------------
example_sce <- runUMAP(example_sce)
head(reducedDim(example_sce, "UMAP"))

## ----plot-reduceddim-4comp-colby-shapeby--------------------------------------
plotReducedDim(example_sce, dimred = "PCA", colour_by = "level1class")

## ----plot-pca-4comp-colby-sizeby-exprs----------------------------------------
plotTSNE(example_sce, colour_by = "Snap25")

## ----plot-pca-default---------------------------------------------------------
plotPCA(example_sce, colour_by="Mog")

## ----plot-pca-4comp-colby-shapeby---------------------------------------------
example_sce <- runPCA(example_sce, ncomponents=20)
plotPCA(example_sce, ncomponents = 4, colour_by = "level1class")

## ---- fig.wide=TRUE-----------------------------------------------------------
ggcells(example_sce, mapping=aes(x=level1class, y=Snap25)) + 
    geom_boxplot() +
    facet_wrap(~tissue)

## -----------------------------------------------------------------------------
ggcells(example_sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=Snap25)) +
    geom_point() +
    stat_density_2d() +
    facet_wrap(~tissue) +
    scale_colour_distiller(direction=1)

## -----------------------------------------------------------------------------
ggcells(example_sce, mapping=aes(x=sizeFactor, y=Actb)) +
    geom_point() +
    geom_smooth()

## -----------------------------------------------------------------------------
colnames(example_sce) <- make.names(colnames(example_sce))
ggfeatures(example_sce, mapping=aes(x=featureType, y=X1772062111_E06)) + 
    geom_violin()

## -----------------------------------------------------------------------------
sessionInfo()

