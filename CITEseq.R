
rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape)
library(corrplot)
library(Matrix)
library(Seurat)
library(SingleR)
library(tidyverse)
library(cowplot)
library(CiteFuse)
library(drake)

 #This workflow starts from a multimodal Seurat object that includes hashtags too, and on which HTO demultiplexing of Seurat and removal of Empty droplets 
# via Droplet Utils has already been run has already been ran. 

#Next steps are to demultiplex data, to call doublets by using hashtags and another computational method, call cell types and do multimodal clustering. 

# Cell typing will be called using singleR on clusters to call non immune cell types and for immune cell types Azimuth will also be sued. We would also do 
# a comparison of the various methods to be added to the markdown report. 

#Finally we will look at the correlation between RNA and protein expression across cell types and overall. 

#-----------------------------
#define functions that are needed and include function library files

#--------------------------------

#============================================================================================================
#main workflow

plan <- drake_plan(
  
 SeuratObj = readRDS("C:/data/10x datasets/To-process/GSE154826-Lung-CITEseq/CITEseq/Citeseq-lung-nonimmune-panel.rds"),
 HTO.annot = 
  
  
  report = target(
    command = {
      rmarkdown::render(knitr_in("analysis-report.Rmd"))
      file_out("analysis-report.html")
      
    }
  )
)


vis_drake_graph(plan)

make(plan)