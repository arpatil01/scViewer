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
library("gridExtra")
library("nebula")
library("dplyr")
library("scater")
library("SingleCellExperiment")

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
                     helpText("Measure the expression of gene of your interest. For ex- PTPRG"), br(), 
                     actionButton("submit1", "Submit", class = "btn-success")),
                     ## primary- color options for buttons
    conditionalPanel(condition="input.tabselected==4", 
                     textInput("gene1",label="First Gene"),
                     helpText("Enter the first gene to measure co-expression. For ex- PTPRG"), br(),
                     textInput("gene2",label="Second Gene"),
                     helpText("Enter the second gene to measure co-expression. For ex- P2RX7"), br(),
                     numericInput("threshold_coexp",label="Threshold", 0), 
                     helpText("Enter the threshold for determining co-expression of genes. \n For ex, Default is 0 i.e., Genes with expression > 0 is said to be expressed."), br(),
                     actionButton("submit_coxp", "Submit", class = "btn-success")),
    conditionalPanel(condition="input.tabselected==5", 
                     textInput("geneGV1",label="Enter Gene Name"), 
                     helpText("Measure the differential expression of gene of your interest. For ex- PTPRG"), br(), 
                     # numericInput("threshold",label="Threshold", 0), 
                     # helpText("Enter the threshold for computing DE of gene- Applied for Fisher's exact test. \n For ex, Default is 0 i.e., Genes with expression > 0 is said to be expressed."), br(),
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
                              and in 'Co-expression Analysis' tab we can visualize the co-expression between two genes of interest at the same time, 
                              and in the 'Differential expression' tab we present differential expression levels between two groups."), br(), br(),
                             h4("Contact: SCD Statistics, Teva Pharmaceuticals", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.4"),
                             print("Abhijeet R. Patil, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br(),
                             print("Gaurav Kumar, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br(),
                             print("Huanyu Zhou, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br(),
                             print("Liling Warren, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"), br(), br()
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

               # h4("UMAP showing cell types", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
               h4("UMAP Plot", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
               print("Below is the UMAP plot showing different cell populations. Each dot is a single-cell and the different colors corresponds to cell types."), br(),
               column(12, align="center", br(), 
                      plotOutput(outputId= 'plot_sum.output', width = "700px", height = "1000px"), #50%
                      
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
                      plotOutput('plot_sum5.Patient.output', width = "1200px", height = "300px"), br(), br(),
                      downloadButton(outputId = "plot_sum5_download.Patient.output",
                                     label = "Download the Co-expression plots for Patient group"), br(),br(),br()),                
               htmlOutput("text3_inputD_1"), br(), 
               
               column(12, align="center", 
                     DT::dataTableOutput(outputId = "plot_sum6_table1.Patient.output", width = "auto", height = "auto"),
                     downloadButton(outputId = "plot_sum6_download.Patient.output",
                                     label = "Download the csv file"), br(), br(), br()),

               htmlOutput("text3_inputD_Patient_Vln"), br(),                
               column(12, align="center", 
                      plotOutput('plot_sum6.Patient_Vln.output', width = "1000px", height = "800px"), br(), br(),
                      downloadButton(outputId = "plot_sum6_download.Patient_Vln.output",
                                     label = "Download the Co-expressing genes Violin plots for Patient group"), br(),br(),br()),                
               
               
               htmlOutput("text3_inputD_2"), br(), 
               column(12, align="center", br(), br(),
                      plotOutput('plot_sum5.Normal.output', width = "1200px", height = "300px"), br(), br(),
                      downloadButton(outputId = "plot_sum5_download.Normal.output",
                                     label = "Download the Co-expression plots for Normal group"), br(),br(),br()),
               
               htmlOutput("text3_inputD_3"), br(), 
               column(12, align="center",
                      DT::dataTableOutput(outputId = "plot_sum6_table1.Normal.output", width = "auto", height = "auto"), #500px 300px
                      downloadButton(outputId = "plot_sum6_download.Normal.output",
                                     label = "Download the csv file"), br(), br(), br()),
               
               htmlOutput("text3_inputD_Normal_Vln"), br(),
               column(12, align="center", 
                      plotOutput('plot_sum6.Normal_Vln.output', width = "1000px", height = "800px"), br(), br(), #800px
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
               ## h5("The Dot plot shows the percent expressed and average expression of the gene in patient and Normal samples across different cell types."), br(),
               column(12, align="center", br(),
                      plotOutput('plot_sum3_DEA.output', width = "800px", height = "400px"), br(),br(), #width = "100%"
                      downloadButton(outputId = "plot_sum3_DEA_download.output",
                                     label = "Download the Dot plot per gene"), br(),br(),br()),
               ## Avg and Perc exp
               htmlOutput("text3a_inputD_2"), br(),
               ## h5("The below Table shows the Average expression and Percent expressed metrics of all the cells irrespective of condition (Case/Control) across different cell types."), br(),
               column(12, align="center", br(),
                      DT::dataTableOutput(outputId = "plot_sum3a_DEA_table1.output",width = "auto", height = "auto"), #1200px 600px
                      downloadButton(outputId = "plot_sum3a_DEA_download.output",
                                     label = "Download the csv file"), br(),br(),br()),
               # ## Fisher's
               # htmlOutput("text4_inputD_2"), br(),
               # column(12, align="center", br(),
               #        DT::dataTableOutput(outputId = "plot_sum4_DEA_table1.output",width = "auto", height = "auto"), #1200px 600px
               #        downloadButton(outputId = "plot_sum4_DEA_download.output",
               #                       label = "Download the csv file"), br(),br(),br()),
               
               ## Wilcoxon
               htmlOutput("text4a_inputD_2"), br(),
               column(12, align="center", br(),
                      DT::dataTableOutput(outputId = "plot_sum4a_DEA_table1.output",width = "auto", height = "auto"), #1200px 600px
                      downloadButton(outputId = "plot_sum4a_DEA_download.output",
                                     label = "Download the csv file"), br(),br(),br()),
               ## Pseudobulk
               htmlOutput("text4_inputD_3"), br(),
               ## h4("Plots showing expression variablity at sample level (Using Pseudo-bulk approach)", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
               column(12, align="center", br(), br(),
                      plotOutput('plot_sum7.output', width = "700px", height = "500px"), br(), br(),
                      downloadButton(outputId = "plot_sum7_download.output",
                                     label = "Download the variability plots"), br(),br(),br())
               ),
      tabPanel(h4("Resources"),
               uiOutput("hlinks")),
      id = "tabselected"
    )
  , width = 10)
)))

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
      # sc_file <- RunUMAP(sc_file, dims = 1:20)
      p1 <- DimPlot(sc_file, reduction = "umap", label=T, label.size=6)
      
      tmp <- data.frame(table(sc_file$cell_type))
      colnames(tmp) <- c("Cell Types", "Number of cells")
      
      p2 <- ggplot(tmp, aes(x=reorder(`Cell Types`, -`Number of cells`), y=`Number of cells`, fill=`Cell Types`))+ 
                   xlab("Cell Types") + ylab("Number of Cells") + #scale_x_discrete(limits = `Number of cells`) +
        geom_bar(stat="identity")+theme_minimal()

      # ggplot(theTable, aes(x=reorder(Position, -table(Position)[Position]))) + geom_bar()
      
      # aes(x=reorder(Position,Position,
      #               function(x)-length(x))))
      
      # p2 <- SCpubr::do_BarPlot(sample=sc_file, 
      #                          group.by = "cell_type",
      #                          # plot.title = "Without split.by - position = fill",
      #                          position = "stack", #colors.use = (hue_pal()(16)),
      #                          flip = FALSE)
      p <- plot_grid(p1, p2,
                     ncol = 1, nrow = 2)
      # p <- CombinePlots(plots = list(p1, p2))
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
       HTML(paste("The below Table shows the Average expression and Percent expressed metrics of all the cells irrespective of condition (Case/Control) across different cell types.  <br>"))
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
                    "The first row shows the gene1, the number of cells, and the percent of cells in which it is expressed, 
                    the second row shows the gene2, the number of cells, and the percent of cells in which it is expressed, 
                    the third row shows the number of cells, and the percent of cells in which both the genes were co-expression.", "<br> <br>"))       }
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
         HTML(paste("<b>","<br> <br> The violin plots show the expression levels of ",gene_query1," and ",gene_query2," genes in the above co-expressing cells (Both) across different cell types in Patient Samples.", "<br>", 
                    "The scatter plot show the correlation of ",gene_query1," and ",gene_query2," genes in the above co-expressing cells","<br>"))         }
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
         HTML(paste("<b>","<br> <br> The violin plots show the expression levels of ",gene_query1," and ",gene_query2," genes in the above co-expressing cells (Both) across different cell types in Normal Samples.", "<br>", 
                    "The scatter plot show the correlation of ",gene_query1," and ",gene_query2," genes in the above co-expressing cells","<br>"))         }
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
         HTML(paste("<b>","Co-expression in Normal samples", "</b>", "<br> <br>",
                    "The Feature plot shows the expression of both genes across different cell types in Normal samples. 
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
         HTML(paste("<b>","The detailed distribution of gene expression among all the cells in Normal samples can be found in below Table.", "</b>", "<br> <br>",
                    "The first row shows the gene1, the number of cells, and the percent of cells in which it is expressed, 
               the second row shows the gene2, the number of cells, and the percent of cells in which it is expressed, 
               the third row shows the number of cells, and the percent of cells in which both the genes were co-expression.", "<br> <br>"))       }
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
        HTML(paste("The Dot plot shows the percent expressed and average expression of the gene in patient and Normal samples across different cell types."))
      }
     }                                   
  })

  ## text3a_input2 <- function()
  ## Avg and Perc Exp 
  text3a_input_2 <- eventReactive(input$submit2,
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
         HTML(paste("The below Table shows the average expression and percent of cells expressing in Case/Control across different cell types.<br>"))
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
         HTML(paste("The below Table shows the differential expression between Case/Control across different cell types using Fisher's exact test. The Fisher's exact test results are based on the threshold's entered.<br>"))
       }
     }                                   
   })
  
  # text4_input2 <- function()
  ## NBGMM
  text4a_input_2 <- eventReactive(input$submit2,
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
       HTML(paste("The below Table shows the differential expression between Case/Control across different cell types using NBGMM tests (Negative Binomial Mixed Models Using Large-Sample Approximation
for Differential Expression Analysis of ScRNA-Seq Data) from nebula package. <br>"))
       # Note: The Wilcoxon rank sum tests does not account for subject level variability.<br>
     }
   }                                   
  })

  ## Pseudobulk-text
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
         HTML(paste("Plots showing average expression at sample level. In each cell type, the average expression of all the individual cells are aggregated at subject-level. <br>"))
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
        p1 <- FeaturePlot(sc_file, gene_query, cols = c("grey","red"), label = T)  #c("grey","red") c("#1BFFFF25", "#2E3192")
        
        # p2 <- SCpubr::do_ViolinPlot(sample = sc_file, features = gene_query)+ xlab("Cell Types")
        p2 <- VlnPlot(sc_file, features = gene_query, pt.size = 0.05)+xlab("Cell Types") #slot = 'scale.data'
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
        # p <- SCpubr::do_DotPlot(sc_file, features = gene_query, group.by=c("cell_type"), legend.type = "colorbar", 
        #                         legend.position = "right", legend.length = 10, rotate_x_axis_labels = 0, scale.by = "size") + 
        #                         xlab('Gene')+ ylab('Cell Types') + 
        #                         ggtitle(paste0("All cells")) + theme(plot.title = element_text(hjust = 0.5))
        
        
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
        
        avg.exp <- DotPlot(sc_file, features=gene_query, group.by=c("cell_type"))
        avg.exp <- avg.exp$data
        avg.exp <- avg.exp[,c(1,2,4)]
        # tmp.case <- round(tmp.case[,c(1,2)],2)
        rownames(avg.exp) <- avg.exp$id
        avg.exp$id <- NULL
        avg.exp <- round(avg.exp, 2)
        colnames(avg.exp) <- c("Average Expression", "Percent Expressed")
        avg.exp <- cbind('Cell Types'=rownames(avg.exp), avg.exp)
        print(avg.exp, row.names = T)
        
        
        
        # avg.exp <- AverageExpression(sc_file, assays = "RNA", gene_query,
        #                              slot = "data", group.by = c("cell_type"))
        # avg.exp <- as.data.frame(t(avg.exp$RNA))
        # avg.exp <- round(avg.exp,2)
        # 
        # avg.exp <- data.frame(names(table(sc_file$cell_type)),
        #                       avg.exp[rownames(avg.exp),])
        # colnames(avg.exp) <- c("Cell Types" , "Average Expression")
        # 
        # ## Subset the cells for percent expressed evaluation
        # perc.exp <- DotPlot(sc_file, features=gene_query, group.by=c("cell_type"))
        # 
        # avg.exp$tmp1 <- perc.exp$data$pct.exp
        # avg.exp$`tmp1` <- round(avg.exp$tmp1,2)
        # colnames(avg.exp)[ncol(avg.exp)] <- paste0("Percent Expressed")
        # 
        # print(avg.exp, row.names = F)
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
        sc_file_AD <- subset(sc_file, subset = grp== paste(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))) ## AD
        sc_file_Normal <- subset(sc_file, subset = grp=="Normal")
        # p1 <- FeaturePlot(sc_file_AD, features=gene_queryGV1,  cols=c("grey","red"), label = T) + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]),"\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
        # p2 <- FeaturePlot(sc_file_Normal, features=gene_queryGV1,  cols=c("grey","red"), label = T)  + ggtitle(paste0("Normal","\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
        # p3 <- VlnPlot(sc_file_AD, features=gene_queryGV1) + xlab('Cell Types')+ ylab('Expression Level') + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]), "\n",gene_queryGV1))  + theme(plot.title = element_text(hjust = 0.5))
        # p4 <- VlnPlot(sc_file_Normal, features=gene_queryGV1) + xlab('Cell Types')+ ylab('Expression Level') + ggtitle(paste0("Normal","\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
        # p <- plot_grid(p1, p2, p3, p4,
        #                ncol = 2, nrow = 2)
        
        # c("grey", "red") c("#1BFFFF25", "#2E3192")
        p1 <- FeaturePlot(object = sc_file, features = gene_queryGV1, split.by = "grp", 
                          keep.scale = "feature", cols = c("grey","red"), label = T) & theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
        p1 <- plot_grid(p1 , ncol = 1, nrow = 1)
        # p2 <- VlnPlot(sc_file, features=gene_queryGV1, split.by = "grp", split.plot = F) + xlab('Cell Types')+ ylab('Expression Level')
        # p <- plot_grid(p1, p2,
        #                ncol = 1, nrow = 2)
        
        p2 <- VlnPlot(sc_file_AD, features=gene_queryGV1, pt.size = 0.05) + xlab('Cell Types')+ ylab('Expression Level') + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]), "\n",gene_queryGV1))  + theme(plot.title = element_text(hjust = 0.5))
        p3 <- VlnPlot(sc_file_Normal, features=gene_queryGV1, pt.size = 0.05) + xlab('Cell Types')+ ylab('Expression Level') + ggtitle(paste0("Normal","\n",gene_queryGV1)) + theme(plot.title = element_text(hjust = 0.5))
        
        
        layout <- rbind(c(1, 1),
                        c(2, 3))
        # layout <- matrix(c(1,1,2,3), 2, 2, byrow = T)
        p <- grid.arrange(p1, p2, p3, layout_matrix=layout)
        
        # p <- grid.arrange(p1,p2,p3, matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
        # p <- plot_grid(p1, p2, p3,
        #                ncol = 2, nrow = 2)
        # layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
        # p1
        # p2
        # p3
        # p <- plot_grid( p1,  p2, p3,
        #   ncol = 2, nrow = 2,
        #   byrow = TRUE
        # )
        # p <- CombinePlots(plots = list(p1, p2, p3))
        
        # p <- ggarrange(p1, p2, p3, #+ rremove("x.text"), 
        #           labels = c("A", "B", "C"),
        #           ncol = 2, nrow = 2)
        
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
        sc_file_AD <- subset(sc_file, subset = grp== paste(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))) ## AD
        sc_file_Normal <- subset(sc_file, subset = grp=="Normal")
        perc.exp_AD <- DotPlot(sc_file_AD, features=gene_queryGV1, group.by=c("cell_type")) + xlab('')+ ylab('Cell Types') + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))) + theme(plot.title = element_text(hjust = 0.5))
        perc.exp_Normal <- DotPlot(sc_file_Normal, features=gene_queryGV1, group.by=c("cell_type")) + xlab('')+ ylab('Cell Types') + ggtitle("Normal") + theme(plot.title = element_text(hjust = 0.5))
        p <- CombinePlots(plots = list(perc.exp_AD, perc.exp_Normal))
        
        # perc.exp_AD <- SCpubr::do_DotPlot(sc_file, features = gene_queryGV1, legend.type = "colorbar",
        #                         legend.position = "right", legend.length = 10,  rotate_x_axis_labels = 0) + xlab('')+ ylab('Cell Types') + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))) + theme(plot.title = element_text(hjust = 0.5))
        # perc.exp_Normal <- SCpubr::do_DotPlot(sc_file, features = gene_queryGV1, legend.type = "colorbar",
        #                         legend.position = "right", legend.length = 10,  rotate_x_axis_labels = 0) + ggtitle(paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))) + theme(plot.title = element_text(hjust = 0.5))
        # p <- CombinePlots(plots = list(perc.exp_AD, perc.exp_Normal))
      }
    }
  }
  )
  
  ## Differential Expression plots-3a (Avg Exp and Perc Exp)
  sum3a_DEA_input <- eventReactive(input$submit2,
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
      thresh = input$threshold
      # thresh=1
      if (input$geneGV1=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
      }
      else
      {

        ## Subset the cells for percent expressed evaluation
        sc_file_AD <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))
        sc_file_Normal <- subset(sc_file, subset = grp=="Normal")
        
        ## Subset the cells for percent expressed evaluation
        tmp.case <- DotPlot(sc_file_AD, features=gene_queryGV1, group.by=c("cell_type"))
        tmp.case <- tmp.case$data
        tmp.case <- tmp.case[,c(1,2,4)]
        # tmp.case <- round(tmp.case[,c(1,2)],2)
        rownames(tmp.case) <- tmp.case$id
        colnames(tmp.case) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", names(table(sc_file_AD$grp))),
                                paste0("Percent Expressed", " in ", names(table(sc_file_AD$grp))),
                                "Cell Types")
        
        tmp.ctl <- DotPlot(sc_file_Normal, features=gene_queryGV1, group.by=c("cell_type"))
        tmp.ctl <- tmp.ctl$data
        tmp.ctl <- tmp.ctl[,c(1,2,4)]
        rownames(tmp.ctl) <- tmp.ctl$id
        colnames(tmp.ctl) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", names(table(sc_file_Normal$grp))),
                                paste0("Percent Expressed", " in ", names(table(sc_file_Normal$grp))),
                                "Cell Types")
        
        df <- merge(tmp.case, tmp.ctl, by="Cell Types")
        
        rownames(df) <- df$`Cell Types`
        df$`Cell Types` <- NULL
        # df[,c(2,3,4,5)]
        df <- round(df, 2)
        # print(df, row.names = F)
        df <- cbind('Cell Types'=rownames(df), df)
        # df$`Cell Types` <- rownames(df)
        print(df, row.names = F)
      }
    }
  }
  )
  
  ## Differential Expression plots-4a (random effects using nebula)
  sum4a_DEA_input <- eventReactive(input$submit2,
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
       # thresh = input$threshold
       # thresh=1
       if (input$geneGV1=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbols", cex = 1.2)
       }
       else
       {
         sce <- as.SingleCellExperiment(sc_file)
         sce.libsize <- librarySizeFactors(sce)

         calc_nebula <- function(df, ex_list, sce.libsize.tmp){
           tryCatch(
             expr = {
               re = nebula(ex_list$count, ex_list$sid, pred=df, offset = sce.libsize.tmp, model = "NBGMM", method = "HL")
               re$summary
               out <- c('P-value'=re$summary$p_cc, 'log2FC'=re$summary$logFC_cc)
             },
             error = function(e){
               message('Looks like groups were not arranged properly!')
               data_g = group_cell(count= ex_list$count,id= ex_list$sid, pred=df)

               is.defined = function(x)is.null(x)

               if (is.defined(data_g))
               {
                 out <- c('P-value'='NA', 'log2FC'='NA') ## If cells already grouped then data_g is NULL useful the gene won't pass filtering we just put NA for that cell type
               }else{
                 re = nebula(data_g$count, data_g$id, pred=data_g$pred, model = "NBGMM", method = "HL")
                 re$summary
                 out <- c('P-value'=re$summary$p_cc, 'log2FC'=re$summary$logFC_cc)
               }
             },
             warning = function(w){
               # message('Caught an warning!')
               print(w)
             },
             finally = {
               # message('All done, quitting.')
             }
           )
         }

         tmp.nbgmm <- list()
         # names(tmp.celltype.stim) <- sort(names(table(sc_file$cell_type)))
         for (i in sort(names(table(sc_file$cell_type))))
         {
           sc_file_tmp <- subset(sc_file, subset = cell_type == i)
           Idents(sc_file_tmp) <- droplevels(Idents(sc_file_tmp))
           sc_file_tmp$cell_type <- droplevels(sc_file_tmp$cell_type)

           cts_mat <- subset(sc_file_tmp, features=c(gene_queryGV1))
           ex_list <- list()

           ## Add count matrix to list
           ex_list$count <- cts_mat@assays$RNA@counts
           # ex_list$count <- cts_mat@assays$SCT@counts

           ## Add sample ids to list
           ex_list$sid <- sc_file_tmp$sample_ids

           ## Add model matrix to list
           pred <- as.data.frame(sc_file_tmp$grp)
           colnames(pred) <- "cc"
           pred$cc = ifelse(pred$cc == "Normal", 0,1)

           df = model.matrix(~cc, data=pred)

           head(sce.libsize[match(colnames(cts_mat),names(sce.libsize))])
           head(colnames(cts_mat))
           length(colnames(cts_mat))
           length(sce.libsize[match(colnames(cts_mat),names(sce.libsize))])
           sce.libsize.tmp <- sce.libsize[match(colnames(cts_mat),names(sce.libsize))]

           tmp.nbgmm[[i]] <- calc_nebula(df, ex_list, sce.libsize.tmp)
         }

         p_val <- vector()
         avg_log2FC <- vector()
         for(i in names(tmp.nbgmm))
         {
           p_val[i] <- tmp.nbgmm[[i]][1]
           avg_log2FC[i] <- tmp.nbgmm[[i]][2]
         }

         nbgmm.df <- data.frame("NBGMM test (P-value)"=unlist(p_val),
                                 "NBGMM test (Avg log2FC)"=unlist(avg_log2FC)
         )
         colnames(nbgmm.df) <- c("NBGMM test (P-value)", "NBGMM test (Avg log2FC)")
         nbgmm.df <- format(nbgmm.df, digits=3)
         df <- cbind('Cell Types' = rownames(nbgmm.df), nbgmm.df)
         print(df, row.names = T)
       }
     }
   }
  )

  ## Pseudobulk based on wilcoxon- average expression values
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

        which_cell <- names(table(sc_file$cell_type))

        ## for 1 gene
        p <- list()
        for (cell.type in which_cell)
        {
          disease_order <- c("Normal", paste0(names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"])))
          t1 <- paste0("^", gene_queryGV1, "$");
          local_cell <- subset(sc_file, subset = cell_type == cell.type)
          gene_var <- c(grep(t1, rownames(local_cell), value = TRUE))

          gene1<- FetchData(local_cell, vars = gene_var)
          colnames(gene1) <- "gene"
          df <- data.frame(local_cell$sample_ids, local_cell$grp, gene1$gene)
          colnames(df) <- c("Samples", "Type", "Gene")

          df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]

          df <- add_column(df, Type = str_extract(str = df$Samples, pattern = "[^_]+"), .after = "Samples")

          p[[cell.type]] <- ggboxplot(df, x = "Type", y = "Gene",
                                      #fill = "Type",
                                      color = "Type",
                                      palette = c("#00AFBB", "#FC4E07"), #, "#E7B800"
                                      # add = c("jitter"), #jitter mean,
                                      order = disease_order,
                                      outlier.shape = NA,
                                      # title = i, round(tmp.val$p.value,5)
                                      ylab = paste("Average Expression"), xlab = "") + ggtitle(paste0(cell.type)) + theme(plot.title = element_text(hjust = 0.5)) + # ggtitle(paste0(cell.type, "\n","P-value = ", format(res$p.value, digits=3))) + #theme(plot.title = element_text(hjust = 0.5))
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + geom_jitter(position=position_jitter(width=.1, height=0))

        }
        plot1 <- CombinePlots(p)
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
        sc_file_AD <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))
        
        p <-FeaturePlot(object = sc_file_AD, features = c(gene_query1, gene_query2), reduction = 'umap',
                        ## pt.size = 1.5, 
                        # pt.size = 0.8,
                        blend = T, cols = c("grey", "red", "blue"), 
                        # blend.threshold = 0.5,
                        label = T)
      }
    }
  }
  )
  
  ## Co-expression analysis
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
        sc_file_Normal <- subset(sc_file, subset = grp=="Normal")
        # cells.use <- WhichCells(sc_file_Normal, expression = gene_query1 > 0 & gene_query2 > 0)
        p <-FeaturePlot(object = sc_file_Normal, features = c(gene_query1, gene_query2), reduction = 'umap',
                        ## pt.size = 1.5, 
                        ## cells.highlight
                        # cells = cells.use,
                        blend = T, cols = c("grey", "red", "blue"), 
                        # pt.size = 0.8,
                        # blend.threshold = 0.5,
                        label = T)
      }
    }
  }
  )
  
  ## Co-expression analysis
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
      thresh_coexp = input$threshold_coexp
      
      if (input$gene1=="" &  input$gene2=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        sc_file <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))
        ## gene co-expression for plotting purposes
        p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                        # pt.size = NULL,
                        blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,

        p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                         # pt.size = 1.5,
                         label = T)
        ## df1
        df1 <- p1$data
        colnames(df1)[4] <- c("GENE1")
        df1$GENE1 <- as.numeric(as.character(df1$GENE1))
        
        p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                         # pt.size = 1.5,
                         label = T)
        ## df2
        df2 <- p2$data
        colnames(df2)[4] <- c("GENE2")
        df2$GENE2 <- as.numeric(as.character(df2$GENE2))
        
        df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
        df[df > thresh_coexp] <- 1 ## Here- Default thresh_coexp=0
        df$Sum <- rowSums(df)
        
        gene1.all <- length(df$GENE1[df$GENE1>0])
        gene2.all <- length(df$GENE2[df$GENE2>0])
        
        genes.both<-table(df$Sum)[names(table(df$Sum))==2]
        names(genes.both) <- NULL
        
        # genes.none<-table(df$Sum)[names(table(df$Sum))==0] ## Here
        # names(genes.none) <- NULL
        
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
        
        # Genes <- c(gene_query1, gene_query2, paste0(gene_query1, " and ", gene_query2, " (Both)")); #, "None"
        Genes <- c(gene_query1, gene_query2, paste0("Both")); #, "None"
        nCells <- c(gene1.all, gene2.all, genes.both);# genes.none
        print(nCells)
        
        #total.cells <- sum(gene1.all, gene2.all, genes.both) # genes.none
        total.cells <- dim(sc_file)[2]
        
        perc.gene1.only <- round(gene1.only/total.cells*100,2)
        perc.gene2.only <- round(gene2.only/total.cells*100,2)
        
        perc.gene1.all <- round(gene1.all/total.cells*100,2)
        perc.gene2.all <- round(gene2.all/total.cells*100,2)
        
        perc.geneboth <- round(genes.both/total.cells*100,2)
        # perc.genenone <- round(genes.none/total.cells*100,2)
        
        # sum(perc.gene1.only, perc.gene2.only, perc.geneboth, perc.genenone)
        # sum(perc.gene1.all, perc.gene2.all, perc.geneboth, perc.genenone)
        
        Percent = c(perc.gene1.all, perc.gene2.all, perc.geneboth) # perc.genenone
        
        df <- data.frame(Genes, nCells, Percent) # , Percent
        print(df, row.names = F)
      }
    }
  }
  )
  ## Co-expression analysis- Violin Plot
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
      thresh_coexp = input$threshold_coexp
      
      if (input$gene1=="" &  input$gene2=="")
      {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
      }
      else
      {
        sc_file <- subset(sc_file, subset = grp==names(table(sc_file$grp)[ names(table(sc_file$grp))!= "Normal"]))
        ## gene co-expression for plotting purposes
        p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                        pt.size = 0.8, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,
        
        p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                         pt.size = 0.8,label = T)
        ## df1
        df1 <- p1$data
        colnames(df1)[4] <- c("GENE1")
        df1$GENE1 <- as.numeric(as.character(df1$GENE1))
        
        p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                         pt.size = 0.8,label = T)
        ## df2
        df2 <- p2$data
        colnames(df2)[4] <- c("GENE2")
        df2$GENE2 <- as.numeric(as.character(df2$GENE2))
        
        df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
        df[df > thresh_coexp] <- 1 
        df$Sum <- rowSums(df)
        
        rownames(df) <- rownames(df1)
        df$GENE1 <- NULL; df$GENE2 <- NULL;
        df$cellnames <- rownames(df)
        
        row_sub <- apply(df, 1, function(row) all(row > 1)) #!=0
        ##Subset as usual
        df <- df[row_sub,]
        
        new_obj <- subset(sc_file, cells = rownames(df))
        
        sp1 <-VlnPlot(new_obj, features = c(gene_query1), pt.size = 0.05) + xlab(paste(names(table(sc_file$grp)), "group")) #, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
        # ggsave(sp1, filename = paste0("coexp_cells_", gene_query1, "_", gene_query1, "_", gene_query2, "_", comp.name, ".pdf"), height = 5, width = 5)
        
        sp2 <-VlnPlot(new_obj, features = c(gene_query2), pt.size = 0.05) + xlab(paste(names(table(sc_file$grp)), "group"))#, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
        # ggsave(sp2, filename = paste0("coexp_cells_", gene_query2, "_", gene_query1, "_", gene_query2, "_", comp.name, ".pdf"), height = 5, width = 5)
        plot_1 <- CombinePlots(plots = list(sp1, sp2))
        
        ## gene expression analysis for individual genes for detailed metrics
        p1 <-FeaturePlot(object = new_obj, features = c(gene_query1), reduction = 'umap',
                         pt.size = 0.8,label = T)
        ## df1
        df1 <- p1$data
        
        p2 <-FeaturePlot(object = new_obj, features = c(gene_query2), reduction = 'umap',
                         pt.size = 0.8,label = T)
        ## df2
        df2 <- p2$data    
        
        plot_list <- list()
        for(i in names(table(sc_file$cell_type)))
        {
          tmpdf1 <- df1[df1$ident==i,]
          Gene1 <- as.numeric(as.character(tmpdf1[,4]))
          
          tmpdf2 <- df2[df2$ident==i,]
          Gene2 <- as.numeric(as.character(tmpdf2[,4]))
          
          ## Combine two data frames
          df <- data.frame(Gene1, Gene2)
          colnames(df) <- c(colnames(df1)[ncol(df1)], colnames(df2)[ncol(df2)])
          
          sp3 <- ggscatter(df, x = colnames(df)[1], y = colnames(df)[2], size = 0.3,
                           add = "reg.line",  # Add regression line
                           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                           conf.int = TRUE # Add confidence interval
          ) + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5))
          # Add correlation coefficient
          plot_list[[i]] <- sp3 + stat_cor(aes(label = ..r.label..),  label.x = 1, 
                                           # method = "kendall", cor.coef.name = "tau",
                                           # method = "pearson"
                                           method = "spearman"
                                           ) #, label.x = 3, label.y = 3 #spearman kendall
        }
        
        plot_2 <- CombinePlots(plots = plot_list)
        p <- plot_grid(plot_1, plot_2, ncol = 1, nrow = 2)
        # p <- CombinePlots(plots = list(sp, sp3))
        
      }
    }
  }
  )

  ## Co-expression analysis
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
       thresh_coexp = input$threshold_coexp
       
       if (input$gene1=="" &  input$gene2=="")
       {
         plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
         text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
       }
       else
       {
         sc_file <- subset(sc_file, subset = grp=="Normal")
         ## gene co-expression for plotting purposes
         p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                         pt.size = 0.8, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,
         
         ## gene expression analysis for individual genes for detailed metrics
         p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                          pt.size = 0.8,label = T)
         ## df1
         df1 <- p1$data
         colnames(df1)[4] <- c("GENE1")
         df1$GENE1 <- as.numeric(as.character(df1$GENE1))
         
         p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                          pt.size = 0.8,label = T)
         ## df2
         df2 <- p2$data
         colnames(df2)[4] <- c("GENE2")
         df2$GENE2 <- as.numeric(as.character(df2$GENE2))
         
         df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
         df[df > thresh_coexp] <- 1 
         df$Sum <- rowSums(df)
         
         gene1.all <- length(df$GENE1[df$GENE1>0])
         gene2.all <- length(df$GENE2[df$GENE2>0])
         
         genes.both<-table(df$Sum)[names(table(df$Sum))==2]
         names(genes.both) <- NULL
         
         # genes.none<-table(df$Sum)[names(table(df$Sum))==0]
         # names(genes.none) <- NULL
         
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
         
         # Genes <- c(gene_query1, gene_query2, paste0(gene_query1, " and ", gene_query2, " (Both)")); #, "None"
         Genes <- c(gene_query1, gene_query2, paste0("Both")); #, "None"
         nCells <- c(gene1.all, gene2.all, genes.both); #genes.none
         print(nCells)
         
         #total.cells <- sum(gene1.all, gene2.all, genes.both) #genes.none
         total.cells <- dim(sc_file)[2]
         
         perc.gene1.only <- round(gene1.only/total.cells*100,2)
         perc.gene2.only <- round(gene2.only/total.cells*100,2)
         
         perc.gene1.all <- round(gene1.all/total.cells*100,2)
         perc.gene2.all <- round(gene2.all/total.cells*100,2)
         
         perc.geneboth <- round(genes.both/total.cells*100,2)
         # perc.genenone <- round(genes.none/total.cells*100,2)
         
         # sum(perc.gene1.only, perc.gene2.only, perc.geneboth, perc.genenone)
         # sum(perc.gene1.all, perc.gene2.all, perc.geneboth, perc.genenone)
         
         Percent = c(perc.gene1.all, perc.gene2.all, perc.geneboth) #perc.genenone
         
         df <- data.frame(Genes, nCells, Percent) # , Percent
         print(df, row.names = F)
         
       }
     }
   }
  )
  
  ## Co-expression analysis Violin plot
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
     thresh_coexp = input$threshold_coexp
     
     if (input$gene1=="" &  input$gene2=="")
     {
       plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
       text(0.5,0.5,"Please enter gene symbol", cex = 1.2)
     }
     else
     {
       sc_file <- subset(sc_file, subset = grp=="Normal")
       ## gene co-expression for plotting purposes
       p <-FeaturePlot(object = sc_file, features = c(gene_query1, gene_query2), reduction = 'umap',
                       pt.size = 0.8, blend = T, cols = c("grey", "red", "blue"), label = T) # blend.threshold = 0.5,
       
       ## gene expression analysis for individual genes for detailed metrics
       p1 <-FeaturePlot(object = sc_file, features = c(gene_query1), reduction = 'umap',
                        pt.size = 0.8,label = T)
       ## df1
       df1 <- p1$data
       colnames(df1)[4] <- c("GENE1")
       df1$GENE1 <- as.numeric(as.character(df1$GENE1))
       
       p2 <-FeaturePlot(object = sc_file, features = c(gene_query2), reduction = 'umap',
                        pt.size = 0.8,label = T)
       ## df2
       df2 <- p2$data
       colnames(df2)[4] <- c("GENE2")
       df2$GENE2 <- as.numeric(as.character(df2$GENE2))
       
       df <- data.frame("GENE1"=df1$GENE1, "GENE2"=df2$GENE2)
       df[df > thresh_coexp] <- 1 
       df$Sum <- rowSums(df)
       
       rownames(df) <- rownames(df1)
       df$GENE1 <- NULL; df$GENE2 <- NULL;
       df$cellnames <- rownames(df)
       
       row_sub <- apply(df, 1, function(row) all(row > 1)) #!=0
       ##Subset as usual
       df <- df[row_sub,]
       
       new_obj <- subset(sc_file, cells = rownames(df))
       
       sp1 <-VlnPlot(new_obj, features = c(gene_query1), pt.size = 0.05) + xlab(paste(names(table(sc_file$grp)), "group")) #, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
       # ggsave(sp1, filename = paste0("coexp_cells_", gene_query1, "_", gene_query1, "_", gene_query2, "_", comp.name, ".pdf"), height = 5, width = 5)
       
       sp2 <-VlnPlot(new_obj, features = c(gene_query2), pt.size = 0.05) + xlab(paste(names(table(sc_file$grp)), "group"))#, c("#e30800", "#f56505", "#006630", "#5b02c7", "#c40080")
       # ggsave(sp2, filename = paste0("coexp_cells_", gene_query2, "_", gene_query1, "_", gene_query2, "_", comp.name, ".pdf"), height = 5, width = 5)
       plot_1 <- CombinePlots(plots = list(sp1, sp2))

       ## gene expression analysis for individual genes for detailed metrics
       p1 <-FeaturePlot(object = new_obj, features = c(gene_query1), reduction = 'umap',
                        pt.size = 0.8,label = T)
       ## df1
       df1 <- p1$data
       
       p2 <-FeaturePlot(object = new_obj, features = c(gene_query2), reduction = 'umap',
                        pt.size = 0.8,label = T)
       ## df2
       df2 <- p2$data    
       
       plot_list <- list()
       for(i in names(table(sc_file$cell_type)))
       {
         tmpdf1 <- df1[df1$ident==i,]
         Gene1 <- as.numeric(as.character(tmpdf1[,4]))
         
         tmpdf2 <- df2[df2$ident==i,]
         Gene2 <- as.numeric(as.character(tmpdf2[,4]))
         
         ## Combine two data frames
         df <- data.frame(Gene1, Gene2)
         colnames(df) <- c(colnames(df1)[ncol(df1)], colnames(df2)[ncol(df2)])
         
         sp3 <- ggscatter(df, x = colnames(df)[1], y = colnames(df)[2], size = 0.3,
                          add = "reg.line",  # Add regression line
                          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE # Add confidence interval
         ) + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5))
         # Add correlation coefficient
         plot_list[[i]] <- sp3 + stat_cor(aes(label = ..r.label..),  label.x = 1, 
                                          # method = "kendall", cor.coef.name = "tau",
                                          # method = "pearson"
                                          method = "spearman"
         ) #, label.x = 3, label.y = 3 #spearman kendall
       }
       
       plot_2 <- CombinePlots(plots = plot_list)
       p <- plot_grid(plot_1, plot_2, ncol = 1, nrow = 2)
       # p <- CombinePlots(plots = list(sp1, sp2, sp3))
       
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
  
  ## Avg exp and Perc exp
  output$text3a_inputD_2 <- renderUI(
  {
    print(text3a_input_2())
  })
  
  # ## Fisher's exact
  # output$text4_inputD_2 <- renderUI(
  # {
  #   print(text4_input_2())
  # })
  
  ## NBGMM
  output$text4a_inputD_2 <- renderUI(
  {
    print(text4a_input_2())
  })

  # Pseudobulk
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
  
  ## Avg and Perc
  output$plot_sum3a_DEA_table1.output <- DT::renderDataTable(
    { 
      datatable(sum3a_DEA_input(), rownames = FALSE)
    }, 
    options = list(autoWidth = TRUE, columnDefs = list(list(width = '100px', targets = "_all")))
  )
  
  # ## Fisher
  # output$plot_sum4_DEA_table1.output <- DT::renderDataTable(
  # { 
  #   datatable(sum4_DEA_input(), rownames = FALSE)
  # }, 
  # options = list(autoWidth = TRUE, columnDefs = list(list(width = '100px', targets = "_all")))
  # )
  
  ## NBGMM
  output$plot_sum4a_DEA_table1.output <- DT::renderDataTable(
  { 
    datatable(sum4a_DEA_input(), rownames = FALSE)
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
    ggsave(file, sum_input(), width = 7, height = 10, units = "in", device = "pdf")
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
  ## Avg and Perc
  output$plot_sum3a_DEA_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste("Avg_Perc_Case_Ctl", '.csv', sep='')
  },
  content= function(file)
  {
    write.csv(sum3a_DEA_input(), file)
  })
  
  # ## Fisher's
  # output$plot_sum4_DEA_download.output <- downloadHandler(
  # filename= function()
  # {
  #   # "Dot_plot.pdf"
  #   paste("Fishers_DEA", '.csv', sep='')
  # },
  # content= function(file)
  # {
  #   write.csv(sum4_DEA_input(), file)
  # })
  
  ## Wilcoxon
  output$plot_sum4a_DEA_download.output <- downloadHandler(
  filename= function()
  {
    # "Dot_plot.pdf"
    paste("NBGMM_DEA", '.csv', sep='')
  },
  content= function(file)
  {
    write.csv(sum4a_DEA_input(), file)
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
    ggsave(file, sum5_input_Patient(), width = 20, height = 5, units = "in", device = "pdf")
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
    ggsave(file, sum5_input_Normal(), width = 20, height = 5, units = "in", device = "pdf")
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
      ggsave(file, sum6_input_Patient_Vln(), width = 10, height = 10,
             units = "in", device = "pdf")
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
    ggsave(file, sum6_input_Normal_Vln(), width = 10, height = 10, 
           units = "in", device = "pdf")
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
    ggsave(file, sum7_input(), width = 10, height = 8, units = "in", device = "pdf")
  })

}))

shinyApp(ui, server)
