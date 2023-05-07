## ----setup, include = FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("reticulate")
module <- py_available() && reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  if (!requireNamespace("devtools"))
#      install.packages("devtools")
#  devtools::install_github("TomKellyGenetics/leiden")

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  install.packages("leiden")

## ---- eval=FALSE, include=FALSE---------------------------------------------------------------------------------------
#  install.packages("leiden",  quiet = TRUE, repos = 1)
#  devtools::install_github("TomKellyGenetics/leiden", ref = "dev")

## ---------------------------------------------------------------------------------------------------------------------
library("leiden")

## ---------------------------------------------------------------------------------------------------------------------
adjacency_matrix <- rbind(cbind(matrix(round(rbinom(400, 1, 0.8)), 20, 20),
                                matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.8)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.2)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.9)), 20, 20)))
str(adjacency_matrix)
dim(adjacency_matrix )

## ---------------------------------------------------------------------------------------------------------------------
library("igraph")
rownames(adjacency_matrix) <- 1:60
colnames(adjacency_matrix) <- 1:60
graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
graph_object

## ---- warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
plot(graph_object, vertex.color = "grey75")

## ---- eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"-------------------------------------------
#  partition <- c(rep(1, 20), rep(2, 20), rep(3, 20))

## ---- eval=module-----------------------------------------------------------------------------------------------------
partition <- leiden(graph_object)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ---- warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ---- eval=module-----------------------------------------------------------------------------------------------------
#run with defaults
  partition <- leiden(graph_object)


#run with ModularityVertexPartition"
  partition <- leiden(graph_object, partition_type = "ModularityVertexPartition")


#run with resolution parameter
  partition <- leiden(graph_object, resolution_parameter = 0.95)

## ---- warning=FALSE, message=FALSE, eval=module-----------------------------------------------------------------------
partition <- leiden(graph_object, resolution_parameter = 0.5)
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ---- warning=FALSE, message=FALSE, eval=module-----------------------------------------------------------------------
partition <- leiden(graph_object, resolution_parameter = 1.8)
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ---- warning=FALSE, message=FALSE, eval=module-----------------------------------------------------------------------
# generate (unweighted) igraph object in R
library("igraph")
adjacency_matrix[adjacency_matrix > 1] <- 1
snn_graph <- graph_from_adjacency_matrix(adjacency_matrix)
partition <- leiden(snn_graph)
table(partition)

## ---- warning=FALSE, message=FALSE, eval=module-----------------------------------------------------------------------
# pass weights to python leidenalg
adjacency_matrix[adjacency_matrix >= 1 ] <- 1
snn_graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted = NULL)
weights <- sample(1:10, sum(adjacency_matrix!=0), replace=TRUE)
partition <- leiden(snn_graph, weights = weights)
table(partition)

## ---- warning=FALSE, message=FALSE, eval=module-----------------------------------------------------------------------
# generate (weighted) igraph object in R
library("igraph")
adjacency_matrix[adjacency_matrix >= 1] <- weights
snn_graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
partition <- leiden(snn_graph)
table(partition)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  library("Seurat")
#  FindClusters(pbmc_small)
#  membership <- leiden(pbmc_small@snn)
#  table(membership)
#  pbmc_small@ident <- as.factor(membership)
#  names(pbmc_small@ident) <- rownames(pbmc_small@meta.data)
#  pbmc_small@meta.data$ident <- as.factor(membership)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  library("Seurat")
#  FindClusters(pbmc_small)
#  membership <- leiden(pbmc_small@graphs$RNA_snn)
#  table(membership)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  FindClusters(pbmc_small, algorithm = "leiden")
#  table(pbmc_small@active.ident)

