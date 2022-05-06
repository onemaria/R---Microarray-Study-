#################################################
# Project Granulocytes/sepsis #
# Maria M
# Dez 2017
#################################################

# There are three groups in the main dataset: A=Sepsis B=Control C=SIRS. 

# The main dataset contains the difference in gene expression between each group was measured. Gene expressions show how an organism 
# responds to different conditions.

# A binary significance value was assigned (1 significant, 0 not significant). 

# Here I looked for genes with significant difference in Sepsis vs SIRS that do not have a significant difference in Control vs SIRS.

# To reveal the specific genes that are upregulated or downregulated due to infection, I selected from siglist a group of genes with 
# the following conditions: a fold change in gene expression between Sepsis vs Control, Sepsis vs SIRS greater than 1 and a significant 
# difference between Control vs SIRS. 

# These genes are the ones that are downregulated or upregulated due to infection

# A second set of genes differentially expressed in Control vs SIRS was selected from siglist. Genes with a fold change in 
# gene expression in Control vs SIRS greater than 2 were selected

# Actions:
# Read siglist
# Extract the genes that have sigAC1=1 & sigBC=0
# Construct Venn Diagrams
# Export Files
# Fold change AB/AC>2 extraction and HEATMAP #
# Fold change AB/AC>1 extraction and HEATMAP #
# Fold change BC>2 extraction and HEATMAP #


library(XLConnect)
library(dplyr)
library(plyr)
library(VennDiagram)
library(gplots)
library(ComplexHeatmap)


############
# Read microarray_data
############

siglist <- read.csv(file="microarray_data.csv", header=TRUE, sep=";")
#total number of genes in microarray_data = 6730

###########
# Extract the genes that have sigAC1=1 & sigBC=0
###########

# select genes that have SigAC=1 (meaning that they are different Sepsis vs SIRS) and BC=0 (meaning they are alike SIRS vs Control conditions)
sig_siglist <-siglist[which(siglist$SigAC==1 & siglist$SigBC ==0),] #2657
write.table(sig_siglist, file = 'sig_siglist.csv', sep=";", dec=".")


###########
# Create Venn Diagrams
###########

#First Diagram: genes with SigAC=1 with genes with SigBC=0
a1 <-filter(siglist, siglist$SigAC==1) #3792
a2 <-filter(siglist, siglist$SigBC==0) #3746

grid.newpage() 
draw.pairwise.venn(nrow(a1), nrow(a2), nrow(sig_siglist),category = c("SigAC=1", "SigBC=0"), 
                   lty = rep("blank",2),fill = c("light blue", "light pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)

#Second Diagram: difAC==1 and p-value<=0.05 and difBC==0 and p-value<=0.05 and
aria11 <-mArray[which(mArray$AdjP_AC <= 0.05),] #3792
aria22 <-mArray[which(mArray$AdjP_BC <= 0.05),] #2984

grid.newpage()
draw.pairwise.venn(nrow(aria11), nrow(aria22), nrow(sig_mArray_p005), category = c("AdjAC p<0.05", "AdjBC p<0.05"), 
                   lty = rep("blank",2), fill = c("light blue", "pink"),alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)


#############
# Export Files
#############
write.csv(mArray, file = 'mArray.csv')
write.csv(sig_mArray, file = 'sig_mArray.csv')
write.csv(sig_mArray_p005, file = 'sig_mArray_p005.csv')


##############################################
# Fold change AB/AC>2 extraction and HEATMAP #
##############################################
list1 <- filter (microarray_data, abs(microarray_data$Diff.of.Group....A.....B.) > 2 & microarray_data$Sig.Index.for.Diff.of.Group....B.....C.== 1)
list2 <- filter (microarray_data, abs(microarray_data$Diff.of.Group....A.....C.) > 2 & microarray_data$Sig.Index.for.Diff.of.Group....B.....C.== 1)
fold_change_2 <- join (list1, list2, by= "GeneSymbol",type = "inner", match = "all")
names(fold_change_2) <- gsub("Diff.of.Group....A.....B.", "Diff.of.Group....A.....B.", names(fold_change_2), ignore.case = TRUE, useBytes = FALSE)
names(fold_change_2) <- gsub("Diff.of.Group....A.....C.", "Diff.of.Group....A.....C.", names(fold_change_2), ignore.case = TRUE, useBytes = FALSE)
names(fold_change_2) <- gsub("Diff.of.Group....B.....C.", "Diff.of.Group....B.....C.", names(fold_change_2), ignore.case = TRUE, useBytes = FALSE)
write.table(fold_change_2, file = 'Fold_change_2_AB_AC.csv', sep=";", dec=".")

heatmap_extraction_2 <- data.frame(GeneSymbol=fold_change_2$GeneSymbol, DiffAB=fold_change_2$DiffAB, DiffAC=fold_change_2$DiffAC, DiffBC=fold_change_2$DiffBC)

my_matrix_1 <-as.matrix(heatmap_extraction_2 [,c(2:4)])
x<- heatmap_extraction_2[,c(1)]
rownames(my_matrix_1)= x
fontsize_1 <- 2

Heatmap(my_matrix_1, row_names_side = "left", row_dend_side = "left", show_column_names = TRUE,
        column_names_side = "top", column_dend_side = "top",
        show_row_names = TRUE, row_names_gp = gpar(cex=fontsize_1), name="Gene Expression")

##############################################
# Fold change AB/AC>1 extraction and HEATMAP #
##############################################

fc_AB <- filter (siglist, abs(siglist$DiffAB) > 1 & siglist$SigBC== 1)
fc_AC <- filter (siglist, abs(siglist$DiffAC) > 1 & siglist$SigBC== 1)
fold_change_AB_AC_1 <- join(fc_AC, fc_AB, by= "GeneSymbol",type = "inner", match = "all")


heatmap_extraction_AB_AC <- data.frame(GeneSymbol=fold_change_AB_AC_1$GeneSymbol, DiffAB=fold_change_AB_AC_1$DiffAB, DiffAC=fold_change_AB_AC_1$DiffAC, DiffBC=fold_change_AB_AC_1$DiffBC)

my_matrix_2 <-as.matrix(heatmap_extraction_AB_AC [,c(2:4)])
x1<- heatmap_extraction_AB_AC[,c(1)]
rownames(my_matrix_2)= x1
fontsize_2 <- 0.25

Heatmap(my_matrix_2, row_names_side = "left", row_dend_side = "left", show_column_names = TRUE,
        column_names_side = "top", column_dend_side = "top",
        show_row_names = TRUE, row_names_gp = gpar(cex=fontsize_2), 
        name="Gene Expression",  width = unit(2.5, "cm"), column_dend_height = unit (1, "cm"), column_names_max_height = unit(3, "cm"))

write.table(fold_change_AB_AC_1, file = 'fold_change_AB_AC_1.csv', sep=";", dec=".")


###########################################
# Fold change BC>2 extraction and HEATMAP #
###########################################

fold_change_BC_2 <- filter (siglist, abs(siglist$DiffBC) > 2)
write.table(Fold_change_BC_2, file = 'Fold_change_BC_2.csv', sep=";", dec=".")

heatmap_extraction_BC <- data.frame(GeneSymbol=fold_change_BC_2$GeneSymbol, DiffAB=fold_change_BC_2$DiffAB, DiffAC=fold_change_BC_2$DiffAC, DiffBC=fold_change_BC_2$DiffBC)

my_matrix <-as.matrix(heatmap_extraction [,c(2:4)])
x2 <- heatmap_extraction_BC[,c(1)]
rownames(my_matrix)= x2
fontsize <- 0.4

Heatmap(my_matrix, row_names_side = "left", row_dend_side = "left", show_column_names = TRUE,
        column_names_side = "top", column_dend_side = "top",
        show_row_names = TRUE, row_names_gp = gpar(cex=fontsize), name="Gene Expression", width = unit(5, "cm"))

