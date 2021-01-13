
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
setwd("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/Rprojects/CITEseq-workflows")
source("CITEseq-functions-seurat-based.R")

#============================================================================================================
#main workflow

plan <- drake_plan(
  
 SeuratObj = readRDS("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/ToTransfer/Citeseq-lung-nonimmune-panel-with-doublet-calls-seurat.rds"),
 Seurat.wnn = Seurat.dimred.wnncluster.CITEseq(SeuratObj),
 SeuratObj.singler = Seurat.Singler(Seurat.wnn),
 SeuratObj.azimuth =Seurat.Azimuth.celltypes(SeuratObj.singler),
    
  report = target(
    command = {
      rmarkdown::render(knitr_in("AnalysisReport.Rmd"))
      file_out("AnalysisReport.html")
      
    }
  )
)


vis_drake_graph(plan)

make(plan)
drake_gc()
