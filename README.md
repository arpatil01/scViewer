# scViewer Application- A single-cell data viewer shiny application

![image](https://user-images.githubusercontent.com/105444693/222185076-9efcc6da-803e-44e8-8fbb-96f72cc0eb22.png)



## Windows users
All the dependencies such as R packages etc are pre-compiled and loaded into the Application bundle
1. Download the scViewer GitHub zip repo and 
2. Unzip the scViewer folder
3. Double click the ```scViewer.bat``` file and the application will be launched

## Linux or Mac users
1. Download the scViewer GitHub zip repo and 
2. Unzip the scViewer folder
3. Install the following packages in your R studio/ R

The package dependencies can also be found here- ```scViewer -> app -> packages.txt```

  ```install.packages("shiny")```
  ```install.packages("shinythemes")```
  ```install.packages("Seurat")```
  ```install.packages("data.table")```
  ```install.packages("DT")```
  ```install.packages("stringr")```
  ```install.packages("ggpubr")```
  ```install.packages("tibble")```

Linux/ Mac users can also redirect by loading the pre-compiled libraries from ```scViewer -> app -> library ``` folder instead of installing them.
For ex, 
library(devtools)
# load package w/o installing using the PATH of pre-installed packages
load_all('/scViewer/app/library/packagename')
 
4. Go to ```scViewer -> app -> shiny``` folder where you can locate the data, www folder, and app.R script to run the shiny app.
5. Run the shiny app using the ```app.R``` script

## Data
The demo dataset can be found in ```scViewer -> app -> data ```

Users can store the processed single-cell RNA_seq data (Seurat .RDS object) in the data folder. The uploaded object will be available in the drop-down menu of Load data tab

## Supplementary Code documentation for processing raw data
The scRNA-seq data processing code can be found in the vignette here ```scViewer -> Supplementary_File.html```

We provide the step-by-step example code documentation for processing the scRNA-seq data. We also provide the metadata format for the scRNA-seq data that is compatible with the app. 
