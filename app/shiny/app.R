## Load required libraries
library("shiny")
library("shinythemes")
library('Seurat')
library('data.table')
library("DT")
library("stringr")
library("ggpubr")
library("cowplot")
library("tibble")
# library('dplyr')
# library("ggplot2")
# library("tidyverse")
# library("gridExtra")
# library("shinyjs")


options(shiny.maxRequestSize = 1000000*1024^2)

rdsfiles <- list.files("./data", pattern = "\\.rds$")

names(rdsfiles) <- gsub(pattern = "\\.rds$", "", rdsfiles)

# studyinfo <- read.csv("./data/studyinfo.csv", check.names = F)

ui <- shinyUI(fluidPage(theme = shinytheme("readable"), pageWithSidebar(
  headerPanel("scViewer- A single-cell RNA sequencing data Viewer"),
  sidebarPanel(
    
    ## conditionalPanel() functions for selected tab
    conditionalPanel(condition="input.tabselected==1",
    ),
    conditionalPanel(condition="input.tabselected==2",
                     selectInput("dataset", "Select a scRNA-seq dataset:", 
                                 # choices = c(" ",rdsfiles)), helpText("Load any pre-processed scRNA-seq data from the above list"), br(), br(),
                                 choices = rdsfiles), 
                     helpText("Load any pre-processed scRNA-seq data from the above list."), 
                     # helpText("For more details about the dataset, Please refer to the Dataset Name column from the Available Dataset Tab."), br(),
                     numericInput("obs", "Number of observations to view from Metadata:",100), helpText("Option to display metadata information")
                     ),
    conditionalPanel(condition="input.tabselected==3", textInput("genes",label="Enter Gene Name"), 
                     helpText("Measure the expression of gene of your interest"), br(), 
                     actionButton("submit1", "Submit", class = "btn-success")),
                     ## primary- color options for buttons
    conditionalPanel(condition="input.tabselected==4", 
                     textInput("gene1",label="First Gene"),
                     helpText("Enter the first gene to measure co-expression"), br(),
                     textInput("gene2",label="Second Gene"),
                     helpText("Enter the first gene to measure co-expression"), br(),
                     actionButton("submit_coxp", "Submit", class = "btn-success")),
    conditionalPanel(condition="input.tabselected==5", 
                     textInput("geneGV1",label="Enter Gene Name"), 
                     helpText("Measure the differential expression of gene of your interest"), br(), 
                     actionButton("submit2", "Submit", class = "btn-success"))
    , width = 2
  ),
  mainPanel(
    # recommend review the syntax for tabsetPanel() & tabPanel() for better understanding 
    # id argument is important in the tabsetPanel()
    # value argument is important in the tabPanle()
    tabsetPanel(
      tabPanel(h4("About"), value = 1, 
               titlePanel(h3("Welcome to scViewer",align="center",style = "font-family: 'times'; color:#081d58; font-si16pt; line-height:1.8", br(), 
                             h4("A single-cell RNA-seq data Viewer Application", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                             # print("Comparison of Single-cell RNA-seq with Bulk RNA-seq (Illustration from 10X genomics)"),br(),
                             # img(src='image2.png', align = "center"), br(), br(),
                             # print("https://www.10xgenomics.com/blog/single-cell-rna-seq-an-introductory-overview-and-tools-for-getting-started"),br(), br(),
                             
                             print("Schematic workflow of the single-cell data Viewer application"),br(),
                             img(src='image1.png', align = "center"), br(), br(),
                             # print("Fresia R, Marangoni P, Burstyn-Cohen T, Sharir A. From Bite to Byte: Dental Structures Resolved at a Single-Cell Resolution. Journal of Dental Research. 2021;100(9):897-905. doi:10.1177/00220345211001848"),br(), br(),
                             
                             ### the rest of your code
                             h4("Input:", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
                             print("A Seurat object in RDS file format"),br(),
                             
                             h4("Output:", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
                             print("It is presented in three tabs, 
                              where in 'Load Data', you will find clustering of your single cell data showing cell types,
                              in 'Expression Plots' tab, we can readily obtain expression profile of your genes of interest across cell types,
                              and in 'Co-expression Analysis' we can visualize the co-expression between two genes of interest at the same time"), br(), br(),
                             h4("Contact: SCD Statistics, Teva Pharmaceuticals", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.4"),
                             print("Abhijeet R. Patil, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br(),
                             print("Gaurav Kumar, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br(),
                             print("Liling Warren, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br(),
                             print("Huanyu Zhou, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br()
               ))),
      
      # tabPanel(h4("Available Datasets"), value = 21,
      #          #########################################################################
      #          h4("List of Available Studies", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
      #          print("Below are the list of single-cell RNA-seq datasets from multiple studies that we have preprocessed and annotated.
      #               Please choose any dataset of your interest from the Dataset Name column below to select in the Load Data Tab."), br(), br(),
      #          DT::dataTableOutput(outputId = "studyinfodata.output"), br(), br()),
      
      tabPanel(h4("Load Data"), value = 2,  
               #########################################################################
               tags$head(tags$style(type="text/css", "#loadmessage {position: fixed;top: 350px;left: 0px;width: 100%;padding: 5px 0px 5px 0px;text-align: center;font-weight: bold;font-size: 100%;color: #000000; z-index: 105;}"),
                         conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                          tags$div("Please wait while the server is loading data...",id="loadmessage"))),
               tags$head(
                 tags$style(HTML(" .shiny-output-error-validation {color: red;}"))
               ), br(), br(),
               print("The fundamental step in the analysis of high-dimensional single-cell data is their visualization in two dimensions. 
               The most widely used nonlinear dimensionality reduction technique is uniform manifold approximation and projection (UMAP). 
               Depending on the study design, the cell types are annotated based on marker genes using manual annotation methods such as- FindMarkers, Differential Expression between one cluster compared to other clusters, or using automatic annotation tools such as scSorter etc. 
               Please refer to the original published study of the loaded dataset to know about cell annotation method applied."),

               h4("UMAP showing cell types", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
               print("Below is the UMAP plot showing different cell populations. Each dot is a single-cell and the different colors corresponds to cell types."), br(),
               column(12, align="center", br(), 
                      plotOutput(outputId= 'plot_sum.output', width = "700px", height = "500px"), #50%
                      
                      downloadButton(outputId = "plot_sum_download.output",
                                     label = "Download the UMAP clustering"), br(), br(),
                      h4("Table showing Metadata information", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(), br(),
                      
                      DT::DTOutput("view"), br(), br())),
      
      tabPanel(h4("Overall Expression"), value = 3, br(), br(),
               htmlOutput("text2_inputD"), br(),
               column(12, align="center", br(), 
                      plotOutput('plot_sum2.output', width = "1000px", height = "500px"), br(),br(), #width = "100%"
                      
                      downloadButton(outputId = "plot_sum2_download.output",
                                     label = "Download the expression plot per gene"), br(),br(),br()
                      ),

               htmlOutput("text2_inputD_1"), br(),
               column(12, align="center", br(),
                      plotOutput('plot_sum3.output', width = "500px", height = "500px"), br(),br(), #width = "100%"
                      downloadButton(outputId = "plot_sum3_download.output",
                                     label = "Download the Dot plot per gene"), br(),br(),br()),
               htmlOutput("text2_inputD_2"), br(),
               column(12, align="center", br(),
                      DT::dataTableOutput(outputId = "plot_sum4_table1.output", width = "500px", height = "400px"),
                      downloadButton(outputId = "plot_sum4_download.output",
                                     label = "Download the csv file"),br())),      
                      
      ## In this tab we will measure the coexpression analysis of a gene
      tabPanel(h4("Co-expression"), value = 4, br(),
               htmlOutput("text3_inputD"), br(),

               column(12, align="center", 
                      plotOutput('plot_sum5.Patient.output', width = "1000px", height = "400px"), br(), br(),
                      downloadButton(outputId = "plot_sum5_download.Patient.output",
                                     label = "Download the Co-expression plots for Patient group"), br(),br(),br()),                
               htmlOutput("text3_inputD_1"), br(), 
               
               column(12, align="center", 
                     DT::dataTableOutput(outputId = "plot_sum6_table1.Patient.output", width = "500px", height = "300px"),
                     downloadButton(outputId = "plot_sum6_download.Patient.output",
                                     label = "Download the csv file"), br(), br(), br()),

               htmlOutput("text3_inputD_Patient_Vln"), br(),                
               column(12, align="center", 
                      plotOutput('plot_sum6.Patient_Vln.output', width = "1000px", height = "400px"), br(), br(),
                      downloadButton(outputId = "plot_sum6_download.Patient_Vln.output",
                                     label = "Download the Co-expressing genes Violin plots for Patient group"), br(),br(),br()),                
               
               
               htmlOutput("text3_inputD_2"), br(), 
               column(12, align="center", br(), br(),
                      plotOutput('plot_sum5.Normal.output', width = "1000px", height = "400px"), br(), br(),
                      downloadButton(outputId = "plot_sum5_download.Normal.output",
                                     label = "Download the Co-expression plots for Normal group"), br(),br(),br()),
               
               htmlOutput("text3_inputD_3"), br(), 
               column(12, align="center",
                      DT::dataTableOutput(outputId = "plot_sum6_table1.Normal.output", width = "500px", height = "300px"),
                      downloadButton(outputId = "plot_sum6_download.Normal.output",
                                     label = "Download the csv file"), br(), br(), br()),
               
               htmlOutput("text3_inputD_Normal_Vln"), br(),
               column(12, align="center", 
                      plotOutput('plot_sum6.Normal_Vln.output', width = "1000px", height = "400px"), br(), br(),
                      downloadButton(outputId = "plot_sum6_download.Normal_Vln.output",
                                     label = "Download the Co-expressing genes Violin plots for Normal group"), br(),br(),br())),
      
      ## In this tab we will validate the gene expression of any or differentially expressed gene using the pseudobulk approach
      tabPanel(h4("Differential Expression"), value = 5, br(), br(),
               htmlOutput("text4_inputD"), br(),
               column(12, align="center", br(), 
                      plotOutput('plot_sum2_DEA.output', width = "800px", height = "800px"), br(),br(),br(), #width = "100%"
                      downloadButton(outputId = "plot_sum2_DEA_download.output", 
                                     label = "Download the expression plot per gene"), br(),br(),br()),
               htmlOutput("text4_inputD_1"), br(),
               # h5("The Dot plot shows the percent expressed and average expression of the gene in patient and normal samples across different cell types."), br(),
               column(12, align="center", br(),
                      plotOutput('plot_sum3_DEA.output', width = "800px", height = "400px"), br(),br(), #width = "100%"
                      downloadButton(outputId = "plot_sum3_DEA_download.output",
                                     label = "Download the Dot plot per gene"), br(),br(),br()),
               htmlOutput("text4_inputD_2"), br(),
               # h5("The below Table shows the Average expression and Percent expressed metrics of all the cells irrespective of condition (Case/Control) across different cell types."), br(),
               column(12, align="center", br(),
                      DT::dataTableOutput(outputId = "plot_sum4_DEA_table1.output",width = "1200px", height = "450px"),
                      downloadButton(outputId = "plot_sum4_DEA_download.output",
                                     label = "Download the csv file"), br(),br(),br()),
               htmlOutput("text4_inputD_3"), br(),
               # h4("Plots showing expression variablity at sample level (Using Pseudo-bulk approach)", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
               column(12, align="center", br(), br(),
                      plotOutput('plot_sum7.output', width = "700px", height = "700px"), br(), br(),
                      downloadButton(outputId = "plot_sum7_download.output",
                                     label = "Download the variability plots"), br(),br(),br())),
      tabPanel(h4("Resources"),
               uiOutput("hlinks")),
      id = "tabselected"
    )
  , width = 10)
)))

onSessionEnded = function(callback) {
  "Registers the given callback to be invoked when the session is closed
(i.e. the connection to the client has been severed). The return value
is a function which unregisters the callback. If multiple callbacks are
registered, the order in which they are invoked is not guaranteed."
  return(.closedCallbacks$register(callback))
}

server <- shinyServer(function(input,output,session)({
  
  session$onSessionEnded(function(){
    stopApp()
  })
  
  # Return the requested dataset
  datasetInput <- reactive(
  {
    df <- readRDS(paste0("./data/", input$dataset))
    return(df)
  })
  
  # Generate a summary of the dataset
  output$summary <- renderPrint(
  {
    dataset <- datasetInput()
    dataset
  })
  
  # Show the first "n" observations
  output$view <- DT::renderDT(
  {
    head(datasetInput(), n = input$obs)
  })
  
  ## Load data- UMAP plot
  sum_input <- reactive(
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      # print("Please select a scRNA-seq dataset for analysis")
      return(invisible())
    }
    else
    {
      sc_file <- rds_file
      
      # umap/tsne
      sc_file <- RunUMAP(sc_file, dims = 1:10)
      p <- DimPlot(sc_file, reduction = "umap", label=T, label.size=6)
    }
  })
  
  # text2_input <- function()
  text2_input <- eventReactive(input$submit1,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      # umap
      gene_query <- toupper(input$genes)
      if (input$genes=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        HTML(paste("Here, we measure the overall expression of a gene in all the cells (Including case and controls) <br>
                   <br> The feature plot and violin plot below the show the expression of the gene in different cell population. <br>"))
      }
    }                                   
  })
  
  # text2_input1 <- function()
  text2_input_1 <- eventReactive(input$submit1,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       # umap
       gene_query <- toupper(input$genes)
       if (input$genes=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         HTML(paste("The Dot plot shows the percent expressed and average expression of the gene in different cell population."))
       }
     }                                   
  })
  
  # text2_input2 <- function()
  text2_input_2 <- eventReactive(input$submit1,
  {
   rds_file <- datasetInput()
   
   if (is.null(rds_file))
   {
     return(invisible())
   }
   else
   { 
     sc_file <- rds_file
     # umap
     gene_query <- toupper(input$genes)
     if (input$genes=="")
     {
       plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
       text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
     }
     else
     {
       HTML(paste("The below Table shows the Average expression and Percent expressed metrics of all the cells irrespective of condition (Case/Control) across different cell types. <br>"))
     }
   }                                   
  })

  
  ## Co-expression analysis
  text3_input <- eventReactive(input$submit_coxp,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query1 <- toupper(input$gene1)
      gene_query2 <- toupper(input$gene2)
      if (input$gene1=="" &  input$gene2=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        HTML(paste("Co-expression analysis is used to explain the functional architecture of genes under different biological conditions.
                     The expression patterns of two genes in a single-cell can be known across different conditions. <br>
                     Here, we can check if two genes are co-expressing in any cell type. <br> <br> <br>", 
                   "<b>","Co-expression in patient samples.", "</b>", "<br> <br>",
                   "The Feature plot shows the expression of both genes across different cell types in patient samples. Each dot in the UMAPs is a cell. 
                      The first feature plot shows the expression of first gene, the second plot shows the expression of second gene, 
                      and the third plot shows the expression levels of both genes across different cell types. 
                      The color levels indicates the expression of gene. The X-axis shows the expression of first gene (Red color), 
                      the Y-axis shows the expression of second gene (Blue color), the cells expressing both genes are shown in gradient of between blue and red (Pink color), 
                      and the cells that do not express either of the genes are marked in grey color.", "<br> <br>"))
      }
    }
  }
  )

  text3_input_1 <- eventReactive(input$submit_coxp,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       
       # umap
       gene_query1 <- toupper(input$gene1)
       gene_query2 <- toupper(input$gene2)
       if (input$gene1=="" &  input$gene2=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         HTML(paste("<b>","The detailed distribution of gene expression among all the cells in patient samples can be found in below Table.", "</b>", "<br> <br>",
                    "The first row shows the gene1 and the number of cells in which it is expressed, the second row shows the gene2 and the number of cells in which it is expressed, 
                      the third row shows the number of cells in which both the genes were co-expression. Finally, none indicates the number of cells among patient samples in which both the genes were not expressed.", "<br> <br>"))       }
     }
   }
  )

  text3_input_Patient_Vln <- eventReactive(input$submit_coxp,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       
       # umap
       gene_query1 <- toupper(input$gene1)
       gene_query2 <- toupper(input$gene2)
       if (input$gene1=="" &  input$gene2=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         HTML(paste("<b>","<br> <br> The violin plot show the expression levels of ",gene_query1," and ",gene_query2," genes in the above co-expressing cells (Both) across different cell types in Patient Samples.", "<br> <br>"))         }
       }
  }
  )  

  text3_input_Normal_Vln <- eventReactive(input$submit_coxp,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       
       # umap
       gene_query1 <- toupper(input$gene1)
       gene_query2 <- toupper(input$gene2)
       if (input$gene1=="" &  input$gene2=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         HTML(paste("<b>","<br> <br> The violin plot show the expression levels of ",gene_query1," and ",gene_query2," genes in the above co-expressing cells (Both) across different cell types in Normal Samples.", "<br> <br>"))}     
       }
  }
  )
      
  text3_input_2 <- eventReactive(input$submit_coxp,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       
       # umap
       gene_query1 <- toupper(input$gene1)
       gene_query2 <- toupper(input$gene2)
       if (input$gene1=="" &  input$gene2=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         HTML(paste("<b>","Co-expression in normal samples", "</b>", "<br> <br>",
                    "The Feature plot shows the expression of both genes across different cell types in normal samples. 
               Each dot in the UMAPs is a cell. The first feature plot shows the expression of first gene, 
               the second plot shows the expression of second gene,
               and the third plot shows the expression levels of both genes across different cell types.
               The color levels indicates the expression of gene. The X-axis shows the expression of first gene (Red color), 
               the Y-axis shows the expression of second gene (Blue color), 
               the cells expressing both genes are shown in gradient of between blue and red (Pink color), 
                     and the cells that do not express either of the genes are marked in grey color.", "<br> <br>"))       
       }
     }
  }
  )
  
  text3_input_3 <- eventReactive(input$submit_coxp,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       
       # umap
       gene_query1 <- toupper(input$gene1)
       gene_query2 <- toupper(input$gene2)
       if (input$gene1=="" &  input$gene2=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         HTML(paste("<b>","The detailed distribution of gene expression among all the cells in normal samples can be found in below Table.", "</b>", "<br> <br>",
                    "The first row shows the gene1 and the number of cells in which it is expressed, 
               the second row shows the gene2 and the number of cells in which it is expressed, 
               the third row shows the number of cells in which both the genes were co-expression. 
               Finally, none indicates the number of cells among normal samples in which both the genes were not expressed.", "<br> <br>"))       }
       }
  }
  )

  text4_input <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
      }
      else
      {
         HTML(paste("Here, we measure the differential expression of a gene in case and control cell population. <br>
          <br> The Feature plots and Violin plots shows the expression of gene between Case and Control cell populations. <br>"))
      }
     }                                   
  })
  
  # text4_input1 <- function()
  text4_input_1 <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
      }
      else
      {
        HTML(paste("The Dot plot shows the percent expressed and average expression of the gene in patient and normal samples across different cell types."))
      }
     }                                   
  })
  
  # text4_input2 <- function()
  text4_input_2 <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
     if (input$geneGV1=="")
     {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
     }
     else
     {
       HTML(paste("The below Table shows the Average expression and Percent expressed metrics of all the cells irrespective of condition (Case/Control) across different cell types. <br>"))
     }
   }                                   
  })

  # text4_input2 <- function()
  text4_input_3 <- eventReactive(input$submit2,
  {
     rds_file <- datasetInput()
     
     if (is.null(rds_file))
     {
       return(invisible())
     }
     else
     { 
       sc_file <- rds_file
       # umap
       gene_queryGV1 <- toupper(input$geneGV1)
       
       if (input$geneGV1=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
       }
       else
       {
         HTML(paste("Plots showing expression variablity at sample level (Using Pseudo-bulk approach) <br>"))
       }
     }                                   
 })
  
    
  ## Expression plots-1
  sum2_input <- eventReactive(input$submit1,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      # umap
      gene_query <- toupper(input$genes)
      if (input$genes=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        p1 <- FeaturePlot(sc_file, gene_query, cols=c("grey","red"), label = T)
        p2 <- VlnPlot(sc_file, features = gene_query)+xlab("Cell Types")
        p <- CombinePlots(plots = list(p1, p2))
      }
    }
  }
  )
  
  ## Expression plots-2
  sum3_input <- eventReactive(input$submit1,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query <- toupper(input$genes)
      if (input$genes=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        perc.exp <- DotPlot(sc_file, features=gene_query, group.by=c("cell_type")) + xlab('Gene')+ ylab('Cell Types') + ggtitle(paste0("All cells")) + theme(plot.title = element_text(hjust = 0.5))
        p <- CombinePlots(plots = list(perc.exp))
      }
    }
  }
  )
  
  ## Expression plots-3
  sum4_input <- eventReactive(input$submit1,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query <- toupper(input$genes)
      if (input$genes=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        avg.exp <- AverageExpression(sc_file, assays = "RNA", gene_query,
                                     slot = "data", group.by = c("cell_type"))
        avg.exp <- as.data.frame(t(avg.exp$RNA))
        avg.exp <- round(avg.exp,2)
        
        avg.exp <- data.frame(names(table(sc_file$cell_type)),
                              avg.exp[rownames(avg.exp),])
        colnames(avg.exp) <- c("Cell Types" , "Average Expression")
        
        ## Subset the cells for percent expressed evaluation
        perc.exp <- DotPlot(sc_file, features=gene_query, group.by=c("cell_type"))
        
        avg.exp$tmp1 <- perc.exp$data$pct.exp
        avg.exp$`tmp1` <- round(avg.exp$tmp1,2)
        colnames(avg.exp)[ncol(avg.exp)] <- paste0("Percent Expressed")
        
        print(avg.exp, row.names = F)
      }
    }
  }
  )
  
  ## Differential Expression plots-1
  sum2_DEA_input <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
            sc_file_AD <- subset(sc_file, subset = grp== paste(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))) ## AD
            sc_file_normal <- subset(sc_file, subset = grp=="normal")
            
            
            p1 <- FeaturePlot(sc_file_AD, features=gene_queryGV1,  cols=c("grey","red"), label = T) + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]),"\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
            p2 <- FeaturePlot(sc_file_normal, features=gene_queryGV1,  cols=c("grey","red"), label = T)  + ggtitle(paste0("Normal","\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
            p3 <- VlnPlot(sc_file_AD, features=gene_queryGV1) + xlab('Cell Types')+ ylab('Expression Level') + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]), "\n",gene_queryGV1))  + theme(plot.title = element_text(hjust = 0.5))
            p4 <- VlnPlot(sc_file_normal, features=gene_queryGV1) + xlab('Cell Types')+ ylab('Expression Level') + ggtitle(paste0("Normal","\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
            p <- plot_grid(p1, p2, p3, p4,
                           ncol = 2, nrow = 2)
      }
    }
  }
  )
  
  ## Differential Expression plots-2
  sum3_DEA_input <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter cell type and gene symbols", cex = 1.2)
      }
      else
      {
        sc_file_AD <- subset(sc_file, subset = grp== paste(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))) ## AD
        sc_file_normal <- subset(sc_file, subset = grp=="normal")
        
        perc.exp_AD <- DotPlot(sc_file_AD, features=gene_queryGV1, group.by=c("cell_type")) + xlab('')+ ylab('Cell Types') + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))) + theme(plot.title = element_text(hjust = 0.5))
        perc.exp_normal <- DotPlot(sc_file_normal, features=gene_queryGV1, group.by=c("cell_type")) + xlab('')+ ylab('Cell Types') + ggtitle("Normal") + theme(plot.title = element_text(hjust = 0.5))
        p <- CombinePlots(plots = list(perc.exp_AD, perc.exp_normal))
      }
    }
  }
  )
  
  ## Differential Expression plots-3
  sum4_DEA_input <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter cell type and gene symbols", cex = 1.2)
      }
      else
      {
        ## Subset the cells for percent expressed evaluation
        sc_file_AD <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))
        sc_file_normal <- subset(sc_file, subset = grp=="normal")
        
        calc_scores <- function(seurat_obj, gene_queryGV1)
        {
          list_AD <- SplitObject(seurat_obj, split.by = "cell_type")
          
          tmp.cellnames.list <- list()
          tmp.cellval.list <- list()
          tmp.cellperc.list <- list()
          tmp.cellcount.list <- list()
          for (p in list_AD) 
          {
            
            tmp.cellnames <- names(table(p$cell_type)[table(p$cell_type)>0])
            tmp.cellnames.list[[length(tmp.cellnames.list)+1]] <- tmp.cellnames
            tmp.cellnames <- unlist(tmp.cellnames.list)
            
            tmp.cellval <- sum(GetAssayData(object = p, slot = "data")[gene_queryGV1,]>0)
            tmp.cellval.list[[length(tmp.cellval.list)+1]] <- tmp.cellval
            tmp.cellval <- unlist(tmp.cellval.list)
            
            tmp.cellcount <- nrow(p@meta.data)
            tmp.cellcount.list[[length(tmp.cellcount.list)+1]] <- tmp.cellcount
            tmp.cellcount <- unlist(tmp.cellcount.list)
            
            tmp.cellperc <- round(sum(GetAssayData(object = p, slot = "data")[gene_queryGV1,]>0)/nrow(p@meta.data)*100,2)
            tmp.cellperc.list[[length(tmp.cellperc.list)+1]] <- tmp.cellperc
            tmp.cellperc <- unlist(tmp.cellperc.list)
            
            tmp.cellval2 <- tmp.cellcount - tmp.cellval
            
            df <- data.frame(tmp.cellnames, tmp.cellval, tmp.cellval2, tmp.cellcount, tmp.cellperc)
            colnames(df) <- c("Cell Types",
                              paste0("# of cells expressing ",gene_queryGV1," in ", names(table(p$grp))),
                              paste0("# of cells not expressing ",gene_queryGV1," in ", names(table(p$grp))),
                              paste0("Total # of cells in ", names(table(p$grp))),
                              paste0("% of cells expressing ",gene_queryGV1," in ", names(table(p$grp))))
          }
          return(df)
        }
        
        calc_fisher_score <- function(df)
        {
          score <- list()
          for(i in 1:nrow(df))
          {
            # print(i)
            tmp1 <- df[i,2]
            tmp2 <- df[i,3]
            tmp3 <- df[i,6]
            tmp4 <- df[i,7]
            tmp.df <- data.frame(c(tmp1,tmp2), c(tmp3, tmp4))
            fish.score <- fisher.test(tmp.df)
            fish.score$p.value <- format(fish.score$p.value, digits=3)
            score[[i]] <- fish.score$p.value
          }
          df$Fisher_Score <- unlist(score)
          colnames(df)[ncol(df)] <- "Fisher's exact test (P-Value)"
          return(df)  
        }
        Dis <- calc_scores(sc_file_AD, gene_queryGV1)
        Normal <- calc_scores(sc_file_normal, gene_queryGV1)
        df <- merge(Dis, Normal, by="Cell Types")
        df <- calc_fisher_score(df)

        print(df, row.names = F)
        
        sc_file$celltype.stim <- paste(Idents(sc_file), sc_file$grp, sep = "_")
        Idents(sc_file) <- "celltype.stim"
        tmp.celltype.stim <- names(table(sc_file$celltype.stim))
        k <- 2
        n <- length(tmp.celltype.stim)
        tmp.celltype.stim <- split(tmp.celltype.stim, rep(1:ceiling(n/k), each=k)[1:n])
        
        tmp.wilcox <- list()
        names(tmp.celltype.stim) <- sort(names(table(sc_file$cell_type)))
        for (i in names(tmp.celltype.stim)) 
        {

          tmp.wilcox[[i]] <- FindMarkers(sc_file, ident.1 = tmp.celltype.stim[[i]][1], ident.2 = tmp.celltype.stim[[i]][2], 
                                  features = gene_queryGV1, 
                                  logfc.threshold=0.0, min.pct = 0, 
                                  min.cells.feature = 0, min.cells.group = 0, 
                                  verbose = FALSE)
        }

        p_val <- vector()
        avg_log2FC <- vector()
        p_val_adj <- vector()
        for(i in names(tmp.wilcox))
        {
          p_val[i] <- tmp.wilcox[[i]][1]
          avg_log2FC[i] <- tmp.wilcox[[i]][2]
          p_val_adj[i] <- tmp.wilcox[[i]][5]
        }
        wilcox.df <- data.frame("Wilcoxon test (P-value)"=unlist(p_val), 
                                "Wilcoxon test (Avg log2FC)"=unlist(avg_log2FC), 
                                "Wilcoxon test (Adj.P-value)"=unlist(p_val_adj))
        colnames(wilcox.df) <- c("Wilcoxon test (P-value)", "Wilcoxon test (Avg log2FC)", "Wilcoxon test (Adj. P-value)")
        rownames(wilcox.df) <- NULL
        wilcox.df <- format(wilcox.df, digits=3)
        df <- cbind(df, wilcox.df)
        print(df, row.names = F)
      }
    }
  }
  )
  
  ## Co-expression analysis
  sum5_input_Patient <- eventReactive(input$submit_coxp,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query1 <- toupper(input$gene1)
      gene_query2 <- toupper(input$gene2)
      if (input$gene1=="" &  input$gene2==""){
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else{
        sc_file_AD <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))
        p <-FeaturePlot(object = sc_file_AD, features = c(gene_query1, gene_query2), reduction = 'umap',
                        pt.size = 1.5, blend = T, cols = c("grey", "red", "blue"), label = T)
      }
    }
  }
  )
  
  sum5_input_Normal <- eventReactive(input$submit_coxp,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query1 <- toupper(input$gene1)
      gene_query2 <- toupper(input$gene2)
      if (input$gene1=="" &  input$gene2==""){
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else{
        sc_file_normal <- subset(sc_file, subset = grp=="normal")
        p <-FeaturePlot(object = sc_file_normal, features = c(gene_query1, gene_query2), reduction = 'umap',
                        pt.size = 1.5, blend = T, cols = c("grey", "red", "blue"), label = T)
      }
    }
  }
  )
  
  sum6_input_Patient <- eventReactive(input$submit_coxp,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query1 <- toupper(input$gene1)
      gene_query2 <- toupper(input$gene2)
      if (input$gene1=="" &  input$gene2=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        sc_file <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))
        ## gene co-expression for plotting purposes
        p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                        pt.size = 1.5, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,

        p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                         pt.size = 1.5,label = T)
        ## df1
        df1 <- p1$data
        colnames(df1)[4] <- c("GENE1")
        df1$GENE1 <- as.numeric(as.character(df1$GENE1))
        
        p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                         pt.size = 1.5,label = T)
        ## df2
        df2 <- p2$data
        colnames(df2)[4] <- c("GENE2")
        df2$GENE2 <- as.numeric(as.character(df2$GENE2))
        
        df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
        df[df > 0] <- 1 
        df$Sum <- rowSums(df)
        
        gene1.all <- length(df$GENE1[df$GENE1>0])
        gene2.all <- length(df$GENE2[df$GENE2>0])
        
        genes.both<-table(df$Sum)[names(table(df$Sum))==2]
        names(genes.both) <- NULL
        
        genes.none<-table(df$Sum)[names(table(df$Sum))==0]
        names(genes.none) <- NULL
        
        gene1.only <- gene1.all-genes.both
        gene2.only <- gene2.all-genes.both
        
        is_empty <- function(x) {
          if (length(x) == 0 & !is.null(x)) {
            x <- 0
          } else {
            x
          }
        }

        genes.both <- is_empty(genes.both)
        
        Genes <- c(gene_query1, gene_query2, paste0(gene_query1, "_", gene_query2, " (Both)"), "None"); 
        nCells <- c(gene1.all, gene2.all, genes.both, genes.none); 
        print(nCells)
        
        df <- data.frame(Genes, nCells) # , Percent
        print(df, row.names = F)
      }
    }
  }
  )

  sum6_input_Patient_Vln <- eventReactive(input$submit_coxp,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query1 <- toupper(input$gene1)
      gene_query2 <- toupper(input$gene2)
      if (input$gene1=="" &  input$gene2=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        sc_file <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"]))
        ## gene co-expression for plotting purposes
        p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                        pt.size = 1.5, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,
        
        p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                         pt.size = 1.5,label = T)
        ## df1
        df1 <- p1$data
        colnames(df1)[4] <- c("GENE1")
        df1$GENE1 <- as.numeric(as.character(df1$GENE1))
        
        p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                         pt.size = 1.5,label = T)
        ## df2
        df2 <- p2$data
        colnames(df2)[4] <- c("GENE2")
        df2$GENE2 <- as.numeric(as.character(df2$GENE2))
        
        df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
        df[df > 0] <- 1 
        df$Sum <- rowSums(df)
        
        rownames(df) <- rownames(df1)
        df$GENE1 <- NULL; df$GENE2 <- NULL;
        df$cellnames <- rownames(df)
        
        row_sub <- apply(df, 1, function(row) all(row > 1)) #!=0
        ##Subset as usual
        df <- df[row_sub,]
        
        new_obj <- subset(sc_file, cells = rownames(df))
        
        sp1 <-VlnPlot(new_obj, features = c(gene_query1), pt.size = 1.5) + xlab(paste(names(table(sc_file$grp)), "group")) #, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
        # ggsave(sp1, filename = paste0("coexp_cells_", list.genes[[i]][1], "_", list.genes[[i]][1], "_", list.genes[[i]][2], "_", comp.name, ".pdf"), height = 5, width = 5)
        
        sp2 <-VlnPlot(new_obj, features = c(gene_query2), pt.size = 1.5) + xlab(paste(names(table(sc_file$grp)), "group"))#, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
        # ggsave(sp2, filename = paste0("coexp_cells_", list.genes[[i]][2], "_", list.genes[[i]][1], "_", list.genes[[i]][2], "_", comp.name, ".pdf"), height = 5, width = 5)
        p <- CombinePlots(plots = list(sp1, sp2))
        
      }
    }
  }
  )

  sum6_input_Normal_Vln <- eventReactive(input$submit_coxp,
  {
   rds_file <- datasetInput()
   
   if (is.null(rds_file))
   {
     return(invisible())
   }
   else
   { 
     sc_file <- rds_file
     
     # umap
     gene_query1 <- toupper(input$gene1)
     gene_query2 <- toupper(input$gene2)
     if (input$gene1=="" &  input$gene2=="")
     {
       plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
       text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
     }
     else
     {
       sc_file <- subset(sc_file, subset = grp=="normal")
       ## gene co-expression for plotting purposes
       p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                       pt.size = 1.5, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,
       
       ## gene expression analysis for individual genes for detailed metrics
       p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                        pt.size = 1.5,label = T)
       ## df1
       df1 <- p1$data
       colnames(df1)[4] <- c("GENE1")
       df1$GENE1 <- as.numeric(as.character(df1$GENE1))
       
       p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                        pt.size = 1.5,label = T)
       ## df2
       df2 <- p2$data
       colnames(df2)[4] <- c("GENE2")
       df2$GENE2 <- as.numeric(as.character(df2$GENE2))
       
       df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
       df[df > 0] <- 1 
       df$Sum <- rowSums(df)
       
       rownames(df) <- rownames(df1)
       df$GENE1 <- NULL; df$GENE2 <- NULL;
       df$cellnames <- rownames(df)
       
       row_sub <- apply(df, 1, function(row) all(row > 1)) #!=0
       ##Subset as usual
       df <- df[row_sub,]
       
       new_obj <- subset(sc_file, cells = rownames(df))
       
       sp1 <-VlnPlot(new_obj, features = c(gene_query1), pt.size = 1.5) + xlab(paste(names(table(sc_file$grp)), "group")) #, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
       # ggsave(sp1, filename = paste0("coexp_cells_", list.genes[[i]][1], "_", list.genes[[i]][1], "_", list.genes[[i]][2], "_", comp.name, ".pdf"), height = 5, width = 5)
       
       sp2 <-VlnPlot(new_obj, features = c(gene_query2), pt.size = 1.5) + xlab(paste(names(table(sc_file$grp)), "group"))#, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
       # ggsave(sp2, filename = paste0("coexp_cells_", list.genes[[i]][2], "_", list.genes[[i]][1], "_", list.genes[[i]][2], "_", comp.name, ".pdf"), height = 5, width = 5)
       p <- CombinePlots(plots = list(sp1, sp2))

     }
   }
  }
  )
  
      
  sum6_input_Normal <- eventReactive(input$submit_coxp,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      
      # umap
      gene_query1 <- toupper(input$gene1)
      gene_query2 <- toupper(input$gene2)
      if (input$gene1=="" &  input$gene2=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        sc_file <- subset(sc_file, subset = grp=="normal")
        ## gene co-expression for plotting purposes
        p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                        pt.size = 1.5, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,
        
        ## gene expression analysis for individual genes for detailed metrics
        p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                         pt.size = 1.5,label = T)
        ## df1
        df1 <- p1$data
        colnames(df1)[4] <- c("GENE1")
        df1$GENE1 <- as.numeric(as.character(df1$GENE1))
        
        p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                         pt.size = 1.5,label = T)
        ## df2
        df2 <- p2$data
        colnames(df2)[4] <- c("GENE2")
        df2$GENE2 <- as.numeric(as.character(df2$GENE2))
        
        df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
        df[df > 0] <- 1 
        df$Sum <- rowSums(df)
        
        gene1.all <- length(df$GENE1[df$GENE1>0])
        gene2.all <- length(df$GENE2[df$GENE2>0])
        
        genes.both<-table(df$Sum)[names(table(df$Sum))==2]
        names(genes.both) <- NULL
        
        genes.none<-table(df$Sum)[names(table(df$Sum))==0]
        names(genes.none) <- NULL
        
        gene1.only <- gene1.all-genes.both
        gene2.only <- gene2.all-genes.both
        
        is_empty <- function(x) {
          if (length(x) == 0 & !is.null(x)) {
            x <- 0
          } else {
            x
          }
        }

        genes.both <- is_empty(genes.both)
        
        Genes <- c(gene_query1, gene_query2, paste0(gene_query1, "_", gene_query2, " (Both)"), "None"); 
        nCells <- c(gene1.all, gene2.all, genes.both, genes.none); 
        print(nCells)
        
        df <- data.frame(Genes, nCells) # , Percent
        print(df, row.names = F)
        
      }
    }
  }
  )
  
  sum7_input <- eventReactive(input$submit2,
  {
    rds_file <- datasetInput()
    
    if (is.null(rds_file))
    {
      return(invisible())
    }
    else
    { 
      sc_file <- rds_file
      # umap
      gene_queryGV1 <- toupper(input$geneGV1)
      
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
      }
      
      else
      {
        # ids <- substr(rownames(sc_file@meta.data), nchar(rownames(sc_file@meta.data))-1,  nchar(rownames(sc_file@meta.data)))
        # sc_file$sample_ids <- paste0(sc_file$grp,ids)
        
        which_cell <- names(table(sc_file$cell_type))
        
        ## for 1 gene
        p <- list()
        for (cell.type in which_cell)
        {
          disease_order <- c("normal", paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "normal"])))
          t1 <- paste0("^", gene_queryGV1, "$"); 
          local_cell <- subset(sc_file, subset = cell_type == cell.type)
          gene_var <- c(grep(t1, rownames(local_cell), value = TRUE))
          
          gene1<- FetchData(local_cell, vars = gene_var)
          colnames(gene1) <- "gene"
          df <- data.frame(local_cell$sample_ids, local_cell$grp, gene1$gene)
          colnames(df) <- c("Samples", "Type", "Gene")
          
          df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
          
          df <- add_column(df, Type = str_extract(str = df$Samples, pattern = "[^_]+"), .after = "Samples")
          res <- wilcox.test(Gene ~ Type, data = df,
                             exact = FALSE)
          # res$p.value
          
          p[[cell.type]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                                     palette = c("#00AFBB", "#FC4E07"), #"#E7B800",
                                     add=c("boxplot", "jitter"),add.params = list(fill="white"),
                                     order = disease_order,
                                     shape = "Type", #size = 0.1,
                                     ylab = gene_var, xlab = FALSE) + ggtitle(paste0(cell.type, "\n","P-value = ", format(res$p.value, digits=3))) + theme(plot.title = element_text(hjust = 0.5))
        }
        plot1 <- CombinePlots(p)
      }
    }
  }
  )
  
  ## URLs
  urldg1 <- a("Trigeminal Ganglion Cell Atlas", href="https://painseq.shinyapps.io/tg-painseq/")
  urldg2 <- a("Pain-seq", href="https://painseq.shinyapps.io/publish/")
  urlhn <- a("Head and Neck", href="https://www.sciencedirect.com/science/article/pii/S0092867417312709?via%3Dihub")
  urlcp <- a("Cross species analysis", href="https://db.cngb.org/nhpca/cross_species")
  urlad <- a("Pre-Frontal Cortex Region", href="https://www.pnas.org/doi/epdf/10.1073/pnas.2008762117")
  urlpd1 <- a("Mid-Brain Region", href="https://academic.oup.com/brain/article/145/3/964/6469020")
  urlpd2 <- a("Pre-Frontal Cortex Region", href="https://www.biorxiv.org/content/10.1101/2022.02.14.480397v1")
  urlcz <- a("CZ CELLXGENE", href="https://cellxgene.cziscience.com/collections")
  urlaz <- a("Azimuth", href="https://azimuth.hubmapconsortium.org/")
  urlbd <- a("Single cell Portal", href="https://singlecell.broadinstitute.org/single_cell")
  urlti <- a("TISCH2", href="http://tisch.comp-genomics.org/")
  urlts <- a("TSGene", href="https://bioinfo.uth.edu/TSGene/?csrt=3941673339932660773")
  urldb <- a("scRNAseqDB", href="https://bioinfo.uth.edu/scrnaseqdb/")
  urlse <- a("CancerSEA", href="http://biocc.hrbmu.edu.cn/CancerSEA/home.jsp")
  output$hlinks <- renderUI(
  {
    tagList(
      h2("Helpful resources for single-cell RNA-seq studies"),
      p("Here we provide some of the additional resources that might be useful"),
      
      h2("Alzhiemer's disease"), 
      p("Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cells and neuroprotective glia in Alzheimer's disease:", urlad),
      
      h2("Parkinson's Disease"),
      p("Single-cell sequencing of human midbrain reveals glial activation and a Parkinson-specific neuronal state:", urlpd1),
      p("Pre-frontal Cortex Region: Single-cell transcriptomic and proteomic analysis of Parkinson's disease Brains:", urlpd2),
      
      h2("Dorsal Root Ganglion"),
      p("Human and Mouse Trigeminal Ganglion Cell Atlas:", urldg1),
      p("Explore Pain Biology at Single-cell Resolution:", urldg2),
      
      h2("Head and Neck"),
      p("Single-cell transcriptomic analysis of primary and metastatic tumor ecosystems in head and neck cancer:", urlhn),
      
      h2("Cross Species Analysis"),
      p("Expore gene expression at single-cell resolution across different species (Human, Mouse, and Monkey):", urlcp),
      
      h2("CZ CELLxGENE"),
      p("single-cell RNA-seq datasets available in cellxgene:", urlcz),
      
      h2("Azimuth"),
      p("Shiny app for reference-based single-cell analysis:", urlaz),
      
      h2("Single cell Portal"),
      p("Single cell research for cancer:", urlbd),
      
      h2("TISCH2"),
      p("Tumor Immune Single-cell Hub 2 (TISCH2) is a scRNA-seq database focusing on tumor microenvironment (TME):", urlti),     
      
      h2("TSGene"),
      p("Tumor Suppressor Gene Database:", urlts),     
      
      h2("scRNAseqDB"),
      p("a database for gene expression profiling in human single cell by RNA-seq:", urldb),  
      
      h2("CancerSEA"),
      p("Cancer Single-cell State Atlas:", urldb),  
    )
  })
  
  ## Render Text
  
  output$text2_inputD <- renderUI(
  {
    print(text2_input())
  })
  
  output$text2_inputD_1 <- renderUI(
  {
    print(text2_input_1())
  })
  
  output$text2_inputD_2 <- renderUI(
  {
    print(text2_input_2())
  })
  
  ############################################
  output$text3_inputD <- renderUI(
  {
    print(text3_input())
  })
  
  output$text3_inputD_1 <- renderUI(
  {
    print(text3_input_1())
  })
  
  output$text3_inputD_Patient_Vln <- renderUI(
  {
    print(text3_input_Patient_Vln())
  })
  
  output$text3_inputD_Normal_Vln <- renderUI(
  {
    print(text3_input_Normal_Vln())
  })
  
  output$text3_inputD_2 <- renderUI(
  {
    print(text3_input_2())
  })
  
  output$text3_inputD_3 <- renderUI(
  {
    print(text3_input_3())
  })
  ############################################
  output$text4_inputD <- renderUI(
  {
    print(text4_input())
  })
  
  output$text4_inputD_1 <- renderUI(
  {
    print(text4_input_1())
  })
  
  output$text4_inputD_2 <- renderUI(
  {
    print(text4_input_2())
  })

  output$text4_inputD_3 <- renderUI(
  {
    print(text4_input_3())
  })
  
  ## Render plots
  output$plot_sum.output <- renderPlot(
  { 
    print(sum_input())
  })
  
  output$plot_sum2.output <- renderPlot(
  { 
    print(sum2_input())
  })
  
  output$plot_sum3.output <- renderPlot(
  { 
    print(sum3_input())
  })
  
  output$plot_sum2_DEA.output <- renderPlot(
  { 
    print(sum2_DEA_input())
  })
  
  output$plot_sum3_DEA.output <- renderPlot(
  { 
    print(sum3_DEA_input())
  })
  
  output$plot_sum5.Patient.output <- renderPlot(
  { 
    print(sum5_input_Patient())
  })
  
  output$plot_sum5.Normal.output <- renderPlot(
  { 
    print(sum5_input_Normal())
  })
  
  output$plot_sum6.Normal_Vln.output <- renderPlot(
  { 
    print(sum6_input_Normal_Vln())
  })
  
  output$plot_sum6.Patient_Vln.output <- renderPlot(
  { 
    print(sum6_input_Patient_Vln())
  })
  
  output$plot_sum7.output <- renderPlot(
  {
    print(sum7_input())
  }
  )
  
  output$plot_sum4_table1.output <- DT::renderDataTable(
  { 
    datatable(sum4_input(), rownames = FALSE)
  }, 
  options = list(autoWidth = TRUE, columnDefs = list(list(width = '200px', targets = "_all")))
  )
  
  output$plot_sum4_DEA_table1.output <- DT::renderDataTable(
  { 
    datatable(sum4_DEA_input(), rownames = FALSE)
  }, 
  options = list(autoWidth = TRUE, columnDefs = list(list(width = '100px', targets = "_all")))
  )
  
  output$plot_sum6_table1.Patient.output <- DT::renderDataTable(
  { 
    datatable(sum6_input_Patient(), rownames = FALSE)
  }, 
  options = list(autoWidth = TRUE, columnDefs = list(list(width = '100px', targets = "_all")))
  )
  
  output$plot_sum6_table1.Normal.output <- DT::renderDataTable(
  { 
    datatable(sum6_input_Normal(), rownames = FALSE)
  }, 
  options = list(autoWidth = TRUE, columnDefs = list(list(width = '100px', targets = "_all")))
  )
  
  ## Downloading plots 
  
  ## Download plot1
  output$plot_sum_download.output <- downloadHandler(
  filename= function()
  {
    # "cellcluster_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  },
  content= function(file)
  {
    ggsave(file, sum_input(), width = 7, height = 7, units = "in", device = "pdf")
  })
  
  ## Download plot2
  output$plot_sum2_download.output <- downloadHandler(
  filename= function()
  {
      paste(input$rdsfile, '.pdf', sep='')
  },
  content= function(file)
  {
    ggsave(file, sum2_input(), width = 10, height = 5, units = "in", device = "pdf")
  })
  
  ## Download plot3
  output$plot_sum3_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum3_input(), width = 7, height = 7, units = "in", device = "pdf")
  })
  
  ## Download table1
  output$plot_sum4_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste("AvgExp_CellTypes", '.csv', sep='')
  },
  content= function(file)
  {
    write.csv(sum4_input(), file)
  })
  
  ## Download plot2_DEA
  output$plot_sum2_DEA_download.output <- downloadHandler(
  filename= function()
  {
    # "Expression_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum2_DEA_input(), width = 12, height = 10, units = "in", device = "pdf")
    
  })
  
  ## Download plot3_DEA
  output$plot_sum3_DEA_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum3_DEA_input(), width = 12, height = 7, units = "in", device = "pdf")
  })
  
  ## Download table1_DEA
  output$plot_sum4_DEA_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste("Gene_CellTypes_DEA", '.csv', sep='')
  },
  content= function(file)
  {
    write.csv(sum4_DEA_input(), file)
  })
  
  ## Download plot4_1
  output$plot_sum5_download.Patient.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum5_input_Patient(), width = 15, height = 5, units = "in", device = "pdf")
  })
  
  ## Download plot4_2
  output$plot_sum5_download.Normal.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum5_input_Normal(), width = 15, height = 5, units = "in", device = "pdf")
  })
  
  ## Download plot_Vln
  output$plot_sum6_download.Patient_Vln.output <- downloadHandler(
    filename= function()
    {
      # "Dot_plot.pdf"
      paste(input$rdsfile, '.pdf', sep='')
    }, 
    content= function(file)
    {
      ggsave(file, sum6_input_Patient_Vln(), width = 8, height = 5, units = "in", device = "pdf")
    })
  
  ## Download plot_Vln
  output$plot_sum6_download.Normal_Vln.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum6_input_Normal_Vln(), width = 8, height = 5, units = "in", device = "pdf")
  })
  
  output$plot_sum6_download.Patient.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste("Coexp_cells_Patients", '.csv', sep='')
  },
  content= function(file)
  {
    write.csv(sum6_input_Patient(), file)
  })
  
  output$plot_sum6_download.Normal.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste("Coexp_cells_Normal", '.csv', sep='')
  },
  content= function(file)
  {
    write.csv(sum6_input_Normal(), file)
  })
  
  ## Download plot5
  output$plot_sum7_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste(input$rdsfile, '.pdf', sep='')
  }, 
  content= function(file)
  {
    ggsave(file, sum7_input(), width = 10, height = 10, units = "in", device = "pdf")
  })

}))

shinyApp(ui, server)