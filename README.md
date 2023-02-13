# scViewer Application- A single-cell data viewer shiny application

![image1](https://user-images.githubusercontent.com/105444693/218345466-95a6b952-34b7-4534-81b7-6b1b1e8c9c86.png)


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
    
4. Go to ```scViewer -> app -> shiny``` folder where you can locate the data, www folder, and app.R script to run the shiny app.
5. Run the shiny app using the ```app.R``` script

## Data
The demo dataset can be found in ```scViewer -> app -> data ```

### Users can store the processed single-cell RNA_seq data (Seurat .RDS object) in the data folder. The uploaded object will be available in the drop-down menu of Load data tab

