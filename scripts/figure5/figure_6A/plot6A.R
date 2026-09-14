library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readxl)


recur_2789.Gene.matrix <- read.csv("../Figure6A.plot.matrix_06_18_25.csv",row.names = 1)
recur_2789_sampleInfo <- read.csv("../Figure6A.plot.info_06_18_25.csv")

#romove samples have no integrations
GeneOnlyOne = apply(recur_2789.Gene.matrix,2,function(col) (sum(col) > 1))
sampleNotZero = apply(recur_2789.Gene.matrix,1,function(row) any(row != 0))

recur_2789.Gene.matrix = recur_2789.Gene.matrix[sampleNotZero,]
recur_2789.Gene.matrix = recur_2789.Gene.matrix[,GeneOnlyOne]
sampleNotZero = apply(recur_2789.Gene.matrix,1,function(row) any(row != 0))

recur_2789.Gene.matrix = recur_2789.Gene.matrix[sampleNotZero,]

recur_2789.Gene.matrix = as.matrix(recur_2789.Gene.matrix)
recur_2789.Gene.matrix.info <- left_join(
  data.frame(Sample = rownames(recur_2789.Gene.matrix)),
  recur_2789_sampleInfo,
  by = c("Sample" = "Sample")
)
n <- length(unique(recur_2789.Gene.matrix.info$patientID))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

patientColor = data.frame(patientID = unique(recur_2789.Gene.matrix.info$patientID),col =colorRampPalette(brewer.pal(12,'Set3'))(n))
temp = as.character(patientColor$col)
names(temp) = as.character(patientColor$patientID)

#remove 075N, 049, 050,074,-21,MIOTO-4497P
#recur_2789.Gene.matrix.info = recur_2789.Gene.matrix.info[-c(32,31,1,2,10,25),]
#row.names.remove <- c('SOP-075N','SOP-049LR3','SOP-050LR3','SOP-074LR1','SOP-021LR1','MIOTO-4497P')
#recur_2789.Gene.matrix = recur_2789.Gene.matrix[!(row.names(recur_2789.Gene.matrix) %in% row.names.remove), ]
recur_2789.Gene.matrix.info$TumorType <- recode_factor(recur_2789.Gene.matrix.info$TumorType, LR = "LR = Local Recurrence", 
                                                       NR = "NR = Node Recurrence","P" = "P = Primary")

pdf('../Figure6A.Gene.heatmap.06.18.25.pdf',width = 20,height = 10)
ha = rowAnnotation(Patient =recur_2789.Gene.matrix.info$patientID,Tumor = recur_2789.Gene.matrix.info$TumorType,
                   HPVtype = recur_2789.Gene.matrix.info$HPVType,
                   col =list(Patient = temp,Tumor = c("LR = Local Recurrence" = "lightblue","NR = Node Recurrence" = "pink"),
                             HPVtype = c("HPV16" = 'orange',"HPV18"= 'blue',"HPV33" = 'green',"None" = 'grey')))
Heatmap(recur_2789.Gene.matrix, 
        name = "Number of HPV integrations", #title of legend
        column_title = "Genes", row_title = "Samples",
        row_title_gp = gpar(fontsize = 20),
        column_title_gp =  gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20), # Text size for row names
        col = circlize::colorRamp2(c(0, 5), c("white", "red")),
        column_dend_height = unit(60, "mm"),
        row_dend_width = unit(60, "mm"),
        show_column_dend = FALSE,
        right_annotation = ha,
        cluster_columns = TRUE,
        clustering_method_columns = 'ward.D',
        cluster_rows = FALSE,
        column_dend_side = c("bottom"))



dev.off()
