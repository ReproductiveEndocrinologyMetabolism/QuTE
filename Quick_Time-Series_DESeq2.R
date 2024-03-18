# Load required packages
library("tximport")
library("readr")
library("tximportData")
library("pasilla")
library("DESeq2")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("cowplot")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("openxlsx")
library("ggvenn")
library("ggsci")
library("scales")
library("clusterProfiler")
show_col(pal_npg("nrc")(6))


# Setting parameters for the analysis
p.value.cutoff = 0.05
p.adj.cutoff = 0.05
p.value.GO = 0.05
q.value.GO = 0.05
gene_count.GO = 3
log2FC.cutoff = 1
convert.Ensembl.GeneSymbol = FALSE
remove.ncRNA = TRUE
do.vst.trans = TRUE
do.log2.trans = FALSE
color.vec = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF")
anno.colours = list(Group = c(PCOS_D0 = color.vec[1], PCOS_D2 = color.vec[2], PCOS_D6 = color.vec[3],
                              Control_D0 = color.vec[4], Control_D2 = color.vec[5], Control_D6 = color.vec[6]))

# Setting colours for plotting
plot.colors = pal_npg("nrc")(6)

# Set which samples (sample ID) to be removed from the count matrix
remove.samples = c() 

# Set genes to be removed or further investigated
remove.genes = c()
PCOS.genes = c()
T2_Diabetes.genes = c()
select.all.genes = c(PCOS.genes, T2_Diabetes.genes)

# Generating the output directory
if (dir.exists(path = paste0("Output/1_DESeq2")) == FALSE) {
  print(paste0("Generating output directory Output/1_DESeq2"))
  dir.create(path = paste0("Output/1_DESeq2"), recursive = TRUE)
  Output.dir = paste0("Output/1_DESeq2/")
} else if (dir.exists(path = paste0("Output/1_DESeq2")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/1_DESeq2/")
} else {
  print("Error with output directory")
}

# Generating the supplementary directory
if (dir.exists(path = paste0("Output/2_Supplementary")) == FALSE) {
  print(paste0("Generating output directory Output/2_Supplementary"))
  dir.create(path = paste0("Output/2_Supplementary"), recursive = TRUE)
  Supp.dir = paste0("Output/2_Supplementary/")
} else if (dir.exists(path = paste0("Output/2_Supplementary")) == TRUE) {
  print("Directory exists")
  Supp.dir = paste0("Output/2_Supplementary/")
} else {
  print("Error with output directory")
}


# Read the data
cts = read.table(file = "Data/countsTabMatrix.txt", sep = "\t", header = TRUE)
rownames(cts) = cts[,1]
cts = cts[,-1]
coldata = openxlsx::read.xlsx(xlsxFile = "Data/Template_experimental_design.xlsx", rowNames = TRUE)
coldata$Group <- factor(coldata$Group)

### Read also ordered coldata for plotting
###coldata.ordered = openxlsx::read.xlsx(xlsxFile = "Data/Template_experimental_design_ordered.xlsx", rowNames = TRUE)
###coldata.ordered$Group <- factor(coldata.ordered$Group)

# Keep only liver samples
cts = cts[, -c(keep.samples)]

# Generate DESeq2 dds object from matrix
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)

# Convert Ensembl to genesymbol, removing ncRNA
if (convert.Ensembl.GeneSymbol == TRUE) {
  dds.ensembl <- rownames(dds)
  dds.symbols <- mapIds(org.Mm.eg.db, keys = dds.ensembl,
                        column = c('SYMBOL'), keytype = 'ENSEMBL')
  dds.symbols <- dds.symbols[!is.na(dds.symbols)]
  dds.symbols <- dds.symbols[match(rownames(dds), names(dds.symbols))]
  rownames(dds) <- dds.symbols
  dds.keep <- !is.na(rownames(dds))
  dds <- dds[dds.keep,]
}

dds

# Run DESeq on the object
dds <- DESeq(dds)

# Extract results PCOS vs. Control each timepoints
res_D0 <- results(dds, contrast=c("Group", "PCOS_D0", "Control_D0"))
res_D2 <- results(dds, contrast=c("Group", "PCOS_D2", "Control_D2"))
res_D6 <- results(dds, contrast=c("Group", "PCOS_D6", "Control_D6"))

# Extract results PCOS over each timepoints
res_PCOS_D0_6 <- results(dds, contrast=c("Group", "PCOS_D2", "PCOS_D0"))
res_PCOS_D2_18 <- results(dds, contrast=c("Group", "PCOS_D6", "PCOS_D2"))
res_PCOS_D0_18 <- results(dds, contrast=c("Group", "PCOS_D6", "PCOS_D0"))

# Extract results Control over each timepoints
res_Ctrl_D0_6 <- results(dds, contrast=c("Group", "Control_D2", "Control_D0"))
res_Ctrl_D2_18 <- results(dds, contrast=c("Group", "Control_D6", "Control_D2"))
res_Ctrl_D0_18 <- results(dds, contrast=c("Group", "Control_D6", "Control_D0"))

# Check summaries of results
summary(res_D0)
summary(res_D2)
summary(res_D6)
summary(res_PCOS_D0_6)
summary(res_PCOS_D2_18)
summary(res_PCOS_D0_18)
summary(res_Ctrl_D0_6)
summary(res_Ctrl_D2_18)
summary(res_Ctrl_D0_18)

# Function to generate, label and output DEG tables
Output.DEG.tables <- function(x, x.name) {
  x.df <- data.frame(x)
  
  # Label the x.df with gene names
  x.ensembl <- rownames(x.df)
  gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys = x.ensembl,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  x.df$GeneName <- gene_symbols
  x.df$EnsemblID <- rownames(x.df)
  x.df$regulation = "Upregulated"
  x.df <- x.df[c("EnsemblID", "GeneName", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "regulation")]
  
  # Label as either up or down regulated
  x.df$regulation[x.df$log2FoldChange < 0] = "Downregulated"
  
  # Output full DE table
  openxlsx::write.xlsx(
    x.df,
    file = paste0(Supp.dir, x.name, "_all_DE_results.xlsx"),
    rowNames = FALSE, overwrite = TRUE
  )
  
  # Filter for p.adj < 0.05, log2FC cutoff > 2
  x.df.05 <- subset(x.df, padj < p.adj.cutoff & (log2FoldChange > log2FC.cutoff | log2FoldChange < -log2FC.cutoff))
  openxlsx::write.xlsx(
    x.df.05,
    file = paste0(Supp.dir, x.name, "_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"_DE_results.xlsx"),
    rowNames = FALSE, overwrite = TRUE
  )
  
  # Remove rows with no gene names
  #x.df.05.no_NA = x.df.05[!is.na(x.df.05$GeneName),]
  #openxlsx::write.xlsx(
  #  x.df.05.no_NA,
  #  file = paste0(Supp.dir, x.name, "_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"NA_removed_DE_results.xlsx"),
  #  rowNames = FALSE, overwrite = TRUE
  #)
  
  # Filter for p.adj < 0.1, log2FC cutoff > 1
  # x.df.1 <- subset(x.df, padj < 0.1 & (log2FoldChange > 2 | log2FoldChange < -2))
  # openxlsx::write.xlsx(x.df.1, file = paste0(Supp.dir,x.name, "_padj0.1_log2FC2_DE_results.xlsx"), rowNames = FALSE)
  
  # Generate vector with significant genes in order
  x.df.05 <- x.df.05[order(x.df.05$log2FoldChange, decreasing = TRUE), ]
  
  return(x.df.05)
}

# Function to generate Volcano plot
Do.VolcanoPlot <- function(x.res, x.name, x.title, volcano.genes = '') {
  
  x.df <- data.frame(x.res)
  
  # Label the x.df with gene names
  x.ensembl <- rownames(x.df)
  gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys = x.ensembl,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  x.df$GeneName <- gene_symbols
  x.df$EnsemblID <- rownames(x.df)
  x.df <- x.df[c("EnsemblID", "GeneName", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  
  # Remove 
  
  # Generate the volcanoplot with gene symbols
  Volcano.plot = EnhancedVolcano(x.df,
                                 lab = x.df$GeneName,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = x.title,
                                 subtitle = '',
                                 ylab = "-Log10 adjusted p-value",
                                 selectLab = volcano.genes,
                                 drawConnectors = TRUE,
                                 arrowheads = FALSE,
                                 legendLabels = '',
                                 legendLabSize = 0,
                                 legendIconSize = 0,
                                 pCutoff = p.adj.cutoff,
                                 FCcutoff = log2FC.cutoff,
                                 pointSize = 2.0,
                                 labSize = 6.0,
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE,
                                 col=c('grey', 'grey', 'grey', 'red3'))
  ggsave2(plot = Volcano.plot, filename = paste0(Output.dir, "Volcano_plot_", x.name, "_significant_DEGs.pdf"), dpi = 700)
  
  
}

# Generate a list to append the DEG tables to
res.list = list("PCOS_Control_D0" = res_D0, "PCOS_Control_D2" = res_D2, "PCOS_Control_D6" = res_D6,
                "PCOS_D0_D2" = res_PCOS_D0_6, "PCOS_D2_D6" = res_PCOS_D2_18, "PCOS_D0_D6" = res_PCOS_D0_18,
                "Ctrl_D0_D2" = res_Ctrl_D0_6, "Ctrl_D2_D6" = res_Ctrl_D2_18, "Ctrl_D0_D6" = res_Ctrl_D0_18)

# Apply Output.DEG.tables and Do.VolcanoPlot to each element in res.list
res.DEG.list = list()
res.titles = c("PCOS D0 vs. Control D0", 
               "PCOS D2 vs. Control D2", 
               "PCOS D6 vs. Control D6",
               "PCOS D0 vs. PCOS D2", 
               "PCOS D2 vs. PCOS D6", 
               "PCOS D0 vs. PCOS D6",
               "Ctrl D0 vs. Ctrl D2", 
               "Ctrl D2 vs. Ctrl D6", 
               "Ctrl D0 vs. Ctrl D6")

for (i in 1:length(res.list)) {
  res.DEG.list[[i]] = Output.DEG.tables(x = res.list[[i]], x.name = names(res.list[i]))
  names(res.DEG.list)[i] = names(res.list[i])
  
  Volcano.plot.DEG = Do.VolcanoPlot(x.res = res.list[[i]], x.name = names(res.list[i]), x.title = res.titles[i])
}

## Extract the genes
res.DEG.genenames = lapply(res.DEG.list, rownames)
names(res.DEG.genenames) = res.titles

# Generate DEG tables corrected for time i.e. removing DEGs in PCOS that are found in Controls
PCOS_D0_D2_corr = res.DEG.list$PCOS_D0_D2[!(res.DEG.list$PCOS_D0_D2$EnsemblID %in% res.DEG.list$Ctrl_D0_D2$EnsemblID), ]
PCOS_D2_D6_corr = res.DEG.list$PCOS_D2_D6[!(res.DEG.list$PCOS_D2_D6$EnsemblID %in% res.DEG.list$Ctrl_D2_D6$EnsemblID), ]
PCOS_D0_D6_corr = res.DEG.list$PCOS_D0_D6[!(res.DEG.list$PCOS_D0_D6$EnsemblID %in% res.DEG.list$Ctrl_D0_D6$EnsemblID), ]

# Extract DEGs that are overlapping
PCOS_D0_D2_overlap = res.DEG.list$PCOS_D0_D2[res.DEG.list$PCOS_D0_D2$EnsemblID %in% res.DEG.list$Ctrl_D0_D2$EnsemblID, ]
PCOS_D2_D6_overlap = res.DEG.list$PCOS_D2_D6[res.DEG.list$PCOS_D2_D6$EnsemblID %in% res.DEG.list$Ctrl_D2_D6$EnsemblID, ]
PCOS_D0_D6_overlap = res.DEG.list$PCOS_D0_D6[res.DEG.list$PCOS_D0_D6$EnsemblID %in% res.DEG.list$Ctrl_D0_D6$EnsemblID, ]

# Attach the corrected and overlapping genes to the result list
res.DEG.list$PCOS_D0_D2_corr = PCOS_D0_D2_corr
res.DEG.list$PCOS_D2_D6_corr = PCOS_D2_D6_corr
res.DEG.list$PCOS_D0_D6_corr = PCOS_D0_D6_corr
res.DEG.list$PCOS_D0_D2_overlap = PCOS_D0_D2_overlap
res.DEG.list$PCOS_D2_D6_overlap = PCOS_D2_D6_overlap
res.DEG.list$PCOS_D0_D6_overlap = PCOS_D0_D6_overlap

# Save the time corrected dataframes
openxlsx::write.xlsx(
  PCOS_D0_D2_corr,
  file = paste0(Supp.dir, "PCOS_D0_D2_padj_", p.adj.cutoff, "_log2FC", log2FC.cutoff,"Time_corrected.xlsx"),
  rowNames = FALSE, overwrite = TRUE
)

openxlsx::write.xlsx(
  PCOS_D2_D6_corr,
  file = paste0(Supp.dir, "PCOS_D2_D6_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"Time_corrected.xlsx"),
  rowNames = FALSE, overwrite = TRUE
)

openxlsx::write.xlsx(
  PCOS_D0_D6_corr,
  file = paste0(Supp.dir, "PCOS_D0_D6_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"Time_corrected.xlsx"),
  rowNames = FALSE, overwrite = TRUE
)

# Update DEG genenames after correction for aging
res.DEG.genenames = lapply(res.DEG.list, rownames)
res.titles = c("PCOS D0 vs. Control D0", 
               "PCOS D2 vs. Control D2", 
               "PCOS D6 vs. Control D6",
               "PCOS D0 vs. PCOS D2", 
               "PCOS D2 vs. PCOS D6", 
               "PCOS D0 vs. PCOS D6",
               "Ctrl D0 vs. Ctrl D2", 
               "Ctrl D2 vs. Ctrl D6", 
               "Ctrl D0 vs. Ctrl D6")
res.titles = c(res.titles, "PCOS D0 vs. PCOS D2 time corr", "PCOS D2 vs. PCOS D6 time corr", "PCOS D0 vs. PCOS D6 time corr",
               "PCOS D0 vs. PCOS D2 overlap", "PCOS D2 vs. PCOS D6 overlap", "PCOS D0 vs. PCOS D6 overlap")

names(res.DEG.genenames) = res.titles



# Generate Venndiagrams of overlapping genes from the three timepoints and corrected for time
## Function to generate the Vennplot and intersect of genes
Do.Vennplot <- function(genelist.x, name.x, Venn.cols = c("#4DBBD5FF", "#E64B35FF", "#00A087FF")) {
  vennplot = ggvenn(genelist.x,
                    fill_color = Venn.cols,
                    stroke_size = 0.5, set_name_size = 5, text_size = 3,
                    show_percentage = FALSE)
  ggsave2(plot = vennplot, filename = paste0(Output.dir, "Venndiagram_", name.x, "_significant_DEGs.pdf"), dpi = 700)
  
  return(vennplot)
}

## Run Vennplot function
Venn.Ctrl_PCOS = Do.Vennplot(genelist.x = res.DEG.genenames[c(1:3)], name.x = "Control_PCOS")
Venn.PCOS = Do.Vennplot(genelist.x = res.DEG.genenames[c(4:6)], name.x = "PCOS")
Venn.Ctrl = Do.Vennplot(genelist.x = res.DEG.genenames[c(7:9)], name.x = "Control")
Venn.PCOS.corr = Do.Vennplot(genelist.x = res.DEG.genenames[c(10:12)], name.x = "PCOS_age_corr")
Venn.PCOS.overlap = Do.Vennplot(genelist.x = res.DEG.genenames[c(13:15)], name.x = "PCOS_age_overlap")
Venn.D0_D2 = Do.Vennplot(genelist.x = res.DEG.genenames[c(7, 4)], name.x = "Control_PCOS_D0_D2")
Venn.D2_D6 = Do.Vennplot(genelist.x = res.DEG.genenames[c(8, 5)], name.x = "Control_PCOS_D2_D6")
Venn_D0_D6 = Do.Vennplot(genelist.x = res.DEG.genenames[c(9, 6)], name.x = "Control_PCOS_D0_D6")

# Merging dataframes to find overlapping and non-overlapping gens
Overlap.DEG.tables <- function(x.df, y.df, z.df, x.label, y.label, z.label, name.x) {
  
  # Overlap by EnsemblID, x.df and y.df
  overlapped.xy = merge(x.df, y.df, 
                        by = "EnsemblID", suffixes = c(paste0("_", x.label), paste0("_", y.label)))
  ensembl.xy = overlapped.xy$EnsemblID
  
  # Overlap by EnsemblID, y.df and z.df
  overlapped.yz = merge(y.df, z.df, 
                        by = "EnsemblID", suffixes = c(paste0("_", y.label), paste0("_", z.label)))
  ensembl.yz = overlapped.yz$EnsemblID
  
  # Overlap by EnsemblID, x.df and z.df
  overlapped.xz = merge(x.df, z.df, 
                        by = "EnsemblID", suffixes = c(paste0("_", x.label), paste0("_", z.label)))
  ensembl.xz = overlapped.xz$EnsemblID
  
  # Overlap by EnsemblID, x.df, y.df and z.df
  ensembl.xyz = intersect(intersect(ensembl.xy, ensembl.xz), ensembl.yz)
  x.df = x.df[x.df$EnsemblID %in% ensembl.xyz,]
  y.df = y.df[y.df$EnsemblID %in% ensembl.xyz,]
  z.df = z.df[z.df$EnsemblID %in% ensembl.xyz,]
  colnames(x.df) <- paste(colnames(x.df), paste0("_", x.label), sep = "")
  colnames(y.df) <- paste(colnames(y.df), paste0("_", y.label), sep = "")
  colnames(z.df) <- paste(colnames(z.df), paste0("_", z.label), sep = "")
  overlapped.xyz = cbind(x.df, y.df, z.df)
  
  # Assign if the gene expression is restored
  overlapped.xy$restored = "No"
  overlapped.xy$restored[overlapped.xy[,9] == "Upregulated" & overlapped.xy[,17] == "Downregulated"] = "Yes"
  overlapped.xy$restored[overlapped.xy[,9] == "Downregulated" & overlapped.xy[,17] == "Upregulated"] = "Yes"
  
  overlapped.yz$restored = "No"
  overlapped.yz$restored[overlapped.yz[,9] == "Upregulated" & overlapped.yz[,17] == "Downregulated"] = "Yes"
  overlapped.yz$restored[overlapped.yz[,9] == "Downregulated" & overlapped.yz[,17] == "Upregulated"] = "Yes"
  
  overlapped.xz$restored = "No"
  overlapped.xz$restored[overlapped.xz[,9] == "Upregulated" & overlapped.xz[,17] == "Downregulated"] = "Yes"
  overlapped.xz$restored[overlapped.xz[,9] == "Downregulated" & overlapped.xz[,17] == "Upregulated"] = "Yes"
  
  # Clean the data matrix
  overlapped.xy = overlapped.xy[,-10]
  names(overlapped.xy)[2] = "GeneName"
  
  overlapped.yz = overlapped.yz[,-10]
  names(overlapped.yz)[2] = "GeneName"
  
  overlapped.xz = overlapped.xz[,-10]
  names(overlapped.xz)[2] = "GeneName"
  
  rownames(overlapped.xyz) = NULL
  overlapped.xyz = overlapped.xyz[,-c(10, 11, 19, 20)]
  names(overlapped.xyz)[2] = "GeneName"
  names(overlapped.xyz)[1] = "EnsemblID"
  
  # Save the overlapped tables in Supplementary 
  openxlsx::write.xlsx(
    overlapped.xy,
    file = paste0(Supp.dir, "Overlapped_", name.x, "_", x.label, "_", y.label, "_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE
  )
  
  openxlsx::write.xlsx(
    overlapped.yz,
    file = paste0(Supp.dir, "Overlapped_", name.x, "_", y.label, "_", z.label, "_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE
  )
  
  openxlsx::write.xlsx(
    overlapped.xz,
    file = paste0(Supp.dir, "Overlapped_", name.x, "_", x.label, "_", z.label, "_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE
  )
  
  openxlsx::write.xlsx(
    overlapped.xyz,
    file = paste0(Supp.dir, "Overlapped_", name.x, "_", x.label, "_", y.label, "_", z.label, "_padj", p.adj.cutoff, "_log2FC", log2FC.cutoff,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE
  )
  
  overlapped.output = list("XY" = overlapped.xy, "YZ" =  overlapped.yz, "XZ" = overlapped.xz, "XYZ" = overlapped.xyz)
  
  # Return the overlapped table
  return(overlapped.output)
  
}

# Run the function and return the overlapped table
res.overlap.all = Overlap.DEG.tables(x.df = res.DEG.list[[1]], y.df = res.DEG.list[[2]], z.df = res.DEG.list[[3]],
                                     x.label = "D0", y.label = "D2", z.label = "D6", name.x = "Control_PCOS")
res.overlap.PCOS = Overlap.DEG.tables(x.df = res.DEG.list[[4]], y.df = res.DEG.list[[5]], z.df = res.DEG.list[[6]],
                                     x.label = "D0_6", y.label = "D2_18", z.label = "D0_18", name.x = "PCOS")
res.overlap.Ctrl = Overlap.DEG.tables(x.df = res.DEG.list[[7]], y.df = res.DEG.list[[8]], z.df = res.DEG.list[[9]],
                                      x.label = "D0_6", y.label = "D2_18", z.label = "D0_18", name.x = "Ctrl")
res.overlap.PCOS_time_corr = Overlap.DEG.tables(x.df = res.DEG.list[[10]], y.df = res.DEG.list[[11]], z.df = res.DEG.list[[12]],
                                              x.label = "D0_6_corr", y.label = "D2_18_corr", z.label = "D0_18_corr", name.x = "PCOS_time_corr")


# If remove.pseudogenes is TRUE, non-coding RNA is removed i.e. keeping only protein coding mRNA
if (remove.ncRNA == TRUE) {
  
  # Remove ncRNA function
  remove.ncRNA.fun <- function(x.df) {
    x.df = x.df[!is.na(x.df$GeneName),]
    return(x.df)
  }
  
  # Run the remove.ncRNA.fun function on all the DEG tables
  for (i in 1:length(res.DEG.list)) {
    res.DEG.list[[i]] = remove.ncRNA.fun(x.df = res.DEG.list[[i]])
    names(res.DEG.list)[i] = names(res.DEG.list[i])
  }
  
}

# New Vennplot after removing ncRNA
res.DEG.genenames.no_ncRNA = lapply(res.DEG.list, rownames)
names(res.DEG.genenames.no_ncRNA) = res.titles

Venn.D0_D2 = Do.Vennplot(genelist.x = res.DEG.genenames.no_ncRNA[c(7, 4)], name.x = "Control_PCOS_D0_D2_no_ncRNA")
Venn.D2_D6 = Do.Vennplot(genelist.x = res.DEG.genenames.no_ncRNA[c(8, 5)], name.x = "Control_PCOS_D2_D6_no_ncRNA")
Venn_D0_D6 = Do.Vennplot(genelist.x = res.DEG.genenames.no_ncRNA[c(9, 6)], name.x = "Control_PCOS_D0_D6_no_ncRNA")

# Perform GO enrichment analysis
Do.GO_enrichment <- function(DEG.df, x.name, number.GOs = 10) {
  
  # Count extract up and down regulated
  up.df = DEG.df[DEG.df$regulation == "Upregulated",]
  down.df = DEG.df[DEG.df$regulation == "Downregulated",]
  
  # Perform GO enrichment analysis
  enrichGO.all = enrichGO(DEG.df$GeneName,
                          OrgDb = org.Mm.eg.db, 
                          ont = "BP",
                          keyType = "SYMBOL",
                          readable = TRUE,
                          pvalueCutoff = p.value.GO,
                          qvalueCutoff = q.value.GO)
  
  enrichGO.all.compare = compareCluster(GeneName~regulation, data=DEG.df, fun="enrichGO",
                                        OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                                        ont = "BP", pAdjustMethod = "BH",
                                        pvalueCutoff  = p.value.GO, qvalueCutoff  = q.value.GO)
  
  enrichGO.up = enrichGO(up.df$GeneName,
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP",
                         keyType = "SYMBOL",
                         readable = TRUE,
                         pvalueCutoff = p.value.GO,
                         qvalueCutoff = q.value.GO)
  
  enrichGO.down = enrichGO(down.df$GeneName,
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP",
                           keyType = "SYMBOL",
                           readable = TRUE,
                           pvalueCutoff = p.value.GO,
                           qvalueCutoff = q.value.GO)
  
  # Extract the result tables
  enrichGO.all.df = enrichGO.all@result
  enrichGO.all.compare.df = enrichGO.all.compare@compareClusterResult
  enrichGO.up.df = enrichGO.up@result
  enrichGO.down.df = enrichGO.down@result
  
  # Save the results tables as .xlsx tables in the supp dir
  openxlsx::write.xlsx(
    enrichGO.all.df,
    file = paste0(Supp.dir, "EnrichGO_BP_All", x.name, "_padj", p.value.GO, "_qval", q.value.GO,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE)
  
  openxlsx::write.xlsx(
    enrichGO.all.compare.df,
    file = paste0(Supp.dir, "EnrichGO_BP_CompareRegulation", x.name, "_padj", p.value.GO, "_qval", q.value.GO,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE)
  
  openxlsx::write.xlsx(
    enrichGO.up.df,
    file = paste0(Supp.dir, "EnrichGO_BP_Upregulated", x.name, "_padj", p.value.GO, "_qval", q.value.GO,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE)
  
  openxlsx::write.xlsx(
    enrichGO.down.df,
    file = paste0(Supp.dir, "EnrichGO_BP_Downregulated", x.name, "_padj", p.value.GO, "_qval", q.value.GO,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE)
  
  # Append to list
  enrichGO.res.list = list("All" = enrichGO.all, "Up" = enrichGO.up, "Down" = enrichGO.down, "Compare" = enrichGO.all.compare)
  
  # Plot the results
  for (i in names(enrichGO.res.list)[1:3]) {
    
    i.dotplot = dotplot(enrichGO.res.list[[i]], showCategory = 20, x = "GeneRatio") + 
      ggtitle(paste0("Dotplot for GSEA BP for ",i, "regulated in ", x.name)) +
      theme_cowplot()
    ggsave2(plot = i.dotplot, filename = paste0(Output.dir, "Dotplot_enrichGO_BP_", x.name, "_", i, "_regulated.pdf"), 
            height = 10, width = 15, dpi = 700)
    
  }
  
  i.dotplot = dotplot(enrichGO.all.compare, showCategory = 10) +
    ggtitle(paste0("Dotplot for GSEA BP for up- and downregulated in ", x.name)) +
    theme_cowplot()
  ggsave2(plot = i.dotplot, filename = paste0(Output.dir, "Dotplot_enrichGO_BP_", x.name, "_UpDown_regulated.pdf"), 
          height = 15, width = 10, dpi = 700)
  
  # Perform -log10 on adjusted p-value
  enrichGO.all.compare.df$Log10_padj = -log10(enrichGO.all.compare.df$p.adjust)
  
  # Make log10_padjust negative for rows where Cluster is "Downregulated"
  enrichGO.all.compare.df$Log10_padj[enrichGO.all.compare.df$Cluster == "Downregulated"] = -abs(enrichGO.all.compare.df$Log10_padj[enrichGO.all.compare.df$Cluster == "Downregulated"])
  
  # Order the dataframe based on log10 p-adjusted
  enrichGO.all.compare.df = enrichGO.all.compare.df[order(enrichGO.all.compare.df$Log10_padj, decreasing = TRUE),]
  
  # Remove GO with less than 3 genes
  enrichGO.all.compare.df = enrichGO.all.compare.df[enrichGO.all.compare.df$Count >= gene_count.GO,]
  enrichGO.all.compare.df = enrichGO.all.compare.df[enrichGO.all.compare.df$p.adjust <= p.value.GO,]
  enrichGO.all.compare.df = enrichGO.all.compare.df[enrichGO.all.compare.df$qvalue <= q.value.GO,]
  
  openxlsx::write.xlsx(
    enrichGO.all.compare.df,
    file = paste0(Supp.dir, "EnrichGO_BP_CompareRegulation_filtered_", x.name, "_padj", p.value.GO, "_qval", q.value.GO,"_table.xlsx"),
    rowNames = FALSE, overwrite = TRUE)
  
  # Extract top and bottom 10
  x.top = head(enrichGO.all.compare.df, number.GOs)
  x.tail = tail(enrichGO.all.compare.df, number.GOs)
  enrichGO.all.compare.top_bottom = rbind(x.top, x.tail)
  
  # Set log10 scale on plot
  log10.min = NULL
  log10.max = NULL
  if (is.null(c(log10.min, log10.max))) {
    log10.min = round(min(enrichGO.all.compare.df$Log10_padj + -10))
    log10.max = round(max(enrichGO.all.compare.df$Log10_padj + 10))
  }
  
  i.dotplot = ggplot(enrichGO.all.compare.top_bottom, aes(x = reorder(Description, Log10_padj), y = Log10_padj, fill = regulation)) + 
    geom_bar(stat = "identity", position = "stack", color = "black") + ylim(log10.min, log10.max) +
    scale_fill_manual(values = c("#0000FF", "#CC0000")) +
    geom_hline(yintercept = 0, color = "black") + coord_flip() + scale_x_discrete(position = "top") + xlab("") + 
    theme_cowplot()
  ggsave2(plot = i.dotplot, filename = paste0(Output.dir, "Barplot_enrichGO_BP_", x.name, "_UpDown_regulated_top", number.GOs, ".pdf"), 
          height = 8, width = 10, dpi = 700)
  
  return(enrichGO.res.list)
  
} 

# # Do GO enrichment for non overlapping genes
GO_PCOS_all_corr_DEGs_D0_6 = Do.GO_enrichment(x.name = "PCOS_all_corr_DEGs_D0_6", DEG.df = res.DEG.list$PCOS_D0_D2_corr)
GO_PCOS_all_corr_DEGs_D2_18 = Do.GO_enrichment(x.name = "PCOS_all_corr_DEGs_D2_18", DEG.df = res.DEG.list$PCOS_D2_D6_corr)
GO_PCOS_all_corr_DEGs_D0_18 = Do.GO_enrichment(x.name = "PCOS_all_corr_DEGs_D0_18", DEG.df = res.DEG.list$PCOS_D0_D6_corr)

# Do GO enrichment for aging overlapping genes
#GO_PCOS_all_overlap_DEGs_D0_6 = Do.GO_enrichment(x.name = "PCOS_all_overlap_DEGs_D0_6", DEG.df = res.DEG.list$PCOS_D0_D2_overlap)
#GO_PCOS_all_overlap_DEGs_D2_18 = Do.GO_enrichment(x.name = "PCOS_all_overlap_DEGs_D2_18", DEG.df = res.DEG.list$PCOS_D2_D6_overlap)
#GO_PCOS_all_overlap_DEGs_D0_18 = Do.GO_enrichment(x.name = "PCOS_all_overlap_DEGs_D0_18", DEG.df = res.DEG.list$PCOS_D0_D6_overlap)

# Perform data transformation on the dds object before heatmap plotting
vst.trans = vst(dds, blind = FALSE)
log2.trans = normTransform(dds)

sd.plot.vst = meanSdPlot(assay(vst.trans))
#ggsave2(plot = sd.plot.vst, filename = paste0(Supp.dir, "MeanSdPlot_VST_transformation.pdf"), dpi = 700)
sd.plot.log2 = meanSdPlot(assay(log2.trans))
#ggsave2(plot = sd.plot.log2, filename = paste0(Supp.dir, "MeanSdPlot_log2_transformation.pdf"), dpi = 700)

# Transform the dds object for downstream plotting
if (do.vst.trans == TRUE && do.log2.trans == FALSE) {
  dds.trans = vst(dds, blind = FALSE)
  trans.x = "vst"
} else if (do.log2.trans == TRUE && do.vst.trans == FALSE) {
  dds.trans = normTransform(dds)
  trans.x = "log2"
}

# Extract transformed groups for downstream plotting of DEGs
res_D0.trans <- dds.trans[, colData(dds.trans)$Group %in% c("PCOS_D0", "Control_D0")]
res_D2.trans <- dds.trans[, colData(dds.trans)$Group %in% c("PCOS_D2", "Control_D2")]
res_D6.trans <- dds.trans[, colData(dds.trans)$Group %in% c("PCOS_D6", "Control_D6")]
res_NIF.trans = dds.trans[, colData(dds.trans)$Group %in% c("PCOS_D0", "PCOS_D2", "PCOS_D6")]
res_Ctrl.trans = dds.trans[, colData(dds.trans)$Group %in% c("Control_D0", "Control_D2", "Control_D6")]

res.trans.list = list("PCOS_Control_D0" = res_D0.trans, "PCOS_Control_D2" = res_D2.trans, "PCOS_Control_D6" = res_D6.trans,
                      "PCOS_D0_D2_D6" = res_NIF.trans, "Control_D0_D2_D6" = res_Ctrl.trans)

# Function to generate heatmaps of transformed objects
Do.Heatmap.fun <- function(x.res, x.name, selected.name, DEG.df = NULL, selected.genes = NULL, top.genes.n = 10, 
                           reorder.col = TRUE, height.x = 10, width.x = 6, selected.break = NULL, set_gaps_cols = c(3,6),
                           convert.ensembl = FALSE) {
  
  # Prepare the transformed object
  #x.df = as.data.frame(assay(x.df)[,order(colnames(assay(x.df)), decreasing = TRUE)])
  x.df = assay(x.res)
  x.df.ensembl = x.df
  rownames(x.df) = mapIds(org.Mm.eg.db, keys = rownames(x.df), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  
  # Remove genes that are NA
  x.df = x.df[!is.na(rownames(x.df)), ]
  
  if (is.null(DEG.df) == TRUE & !is.null(selected.genes) == TRUE) {
    
    # If true, ensembl is converted to gene names
    if (convert.ensembl == TRUE) {
      selected.genes = mapIds(org.Mm.eg.db, keys = selected.genes, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
      selected.genes = selected.genes[!is.na(selected.genes)]
      selected.genes.index = match(selected.genes, rownames(x.df))
    } else if (convert.ensembl == FALSE) {
      selected.genes = selected.genes[!is.na(selected.genes)]
      selected.genes.index = match(selected.genes, rownames(x.df))
    }
    
    
  } else if (!is.null(DEG.df) == TRUE & is.null(selected.genes) == TRUE) {
    DEG.df = DEG.df[!is.na(DEG.df$GeneName),]
  }
  
  
  if (reorder.col == TRUE) {
    x.df = as.data.frame(x.df[,order(colnames(assay(x.res)), decreasing = TRUE)])
  }
  
  # Do heatmap on top and all genes
  if (is.null(selected.genes) == TRUE) {
    
    # Count how many up and downregulated
    up.n = nrow(DEG.df[DEG.df$regulation == "Upregulated",])
    down.n = nrow(DEG.df[DEG.df$regulation == "Downregulated",])
    
    # Extract top n genes
    if (up.n >= top.genes.n && down.n >= top.genes.n) {
      top.genes = DEG.df[c(1:top.genes.n, (nrow(DEG.df)-top.genes.n+1):(nrow(DEG.df))),]
      break.n = top.genes.n
    } else if (up.n >= top.genes.n && down.n < top.genes.n) {
      top.genes = DEG.df[c(1:top.genes.n, (nrow(DEG.df)-down.n+1):(nrow(DEG.df))),]
      break.n = top.genes.n
    } else if (up.n < top.genes.n && down.n >= top.genes.n) {
      top.genes = DEG.df[c(1:up.n, (nrow(DEG.df)-top.genes.n+1):(nrow(DEG.df))),]
      break.n = up.n
    } else if (up.n < top.genes.n && down.n < top.genes.n) {
      top.genes = DEG.df
      break.n = up.n
    }
    
    # Extract the index of top and all genes in the transformed object
    top.genes.index = match(top.genes$GeneName, rownames(x.df))
    top.genes.index = top.genes.index[!is.na(top.genes.index)]
    
    all.genes.index = match(DEG.df$GeneName, rownames(x.df))
    all.genes.index = all.genes.index[!is.na(all.genes.index)]
    
    # Generate the heatmap of top DEGs
    heatmap = pheatmap(x.df[top.genes.index,], cluster_rows=FALSE, 
                       show_rownames=TRUE, show_colnames = TRUE,
                       cluster_cols=FALSE, annotation_col = coldata,
                       annotation_colors = anno.colours, 
                       gaps_col = set_gaps_cols, gaps_row = break.n, scale = 'row')
    ggsave2(plot = heatmap, filename = paste0(Output.dir, "Heatmap_", x.name, "_top_", top.genes.n, ".pdf"), 
            height = height.x, width = width.x, dpi = 700)
    
    # Generate the heatmap of top DEGs
    heatmap = pheatmap(x.df[top.genes.index,], cluster_rows=TRUE,
                       clustering_distance_rows = "euclidean",
                       gaps_col = set_gaps_cols, gaps_row = break.n,
                       show_rownames=TRUE, show_colnames = TRUE,
                       annotation_colors = anno.colours,
                       cluster_cols=FALSE, annotation_col = coldata, scale = 'row')
    ggsave2(plot = heatmap, filename = paste0(Output.dir, "Heatmap_", x.name, "_top_clustered", top.genes.n, ".pdf"), 
            height = height.x, width = width.x, dpi = 700)
    
    # Generate the heatmaps of all DEGs
    x.all.df = x.df[all.genes.index,]
    x.all.df = na.omit(x.all.df)
    
    heatmap = pheatmap(x.all.df, cluster_rows=FALSE, 
                       show_rownames=FALSE,show_colnames = TRUE,
                       cluster_cols=FALSE, annotation_col = coldata,
                       annotation_colors = anno.colours, gaps_col = set_gaps_cols, scale = 'row')
    ggsave2(plot = heatmap, filename = paste0(Output.dir, "Heatmap_", x.name, "_all_DEGs.pdf"), 
            height = height.x, width = width.x, dpi = 700)
    
  } else if (is.null(selected.genes) == FALSE) {
    
    # Generate the heatmap of the selected genes
    heatmap = pheatmap(x.df[selected.genes.index,], cluster_rows=FALSE, 
                       show_rownames=TRUE,show_colnames = TRUE,
                       cluster_cols=FALSE, annotation_col = coldata,
                       annotation_colors = anno.colours, 
                       gaps_col = set_gaps_cols, 
                       gaps_row = selected.break,
                       scale = 'row')
    ggsave2(plot = heatmap, filename = paste0(Output.dir, "Heatmap_", x.name, "_", selected.name, "_selected_genes.pdf"), 
            height = height.x, width = width.x, dpi = 700)
    
    heatmap = pheatmap(x.df[selected.genes.index,], cluster_rows=TRUE, 
                       show_rownames=TRUE,show_colnames = TRUE,
                       cluster_cols=FALSE, annotation_col = coldata,
                       annotation_colors = anno.colours,
                       scale = 'row')
    ggsave2(plot = heatmap, filename = paste0(Output.dir, "Heatmap_", x.name, "_", selected.name, "_selected_genes_clustered.pdf"), 
            height = height.x, width = width.x, dpi = 700)
    
  }
  
  
}

# Do heatmaps two groups, Ctrl vs. PCOS
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[1]], x.name = names(res.trans.list[1]), DEG.df = res.DEG.list$PCOS_Control_D0, 
                             top.genes.n = 25, reorder.col = TRUE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[2]], x.name = names(res.trans.list[2]), DEG.df = res.DEG.list$PCOS_Control_D2, 
                             top.genes.n = 25, reorder.col = TRUE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[3]], x.name = names(res.trans.list[3]), DEG.df = res.DEG.list$PCOS_Control_D6, 
                             top.genes.n = 25, reorder.col = TRUE, set_gaps_cols = 2)

# Do heatmaps, three groups over time, PCOS
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_DEGs_D0_6", DEG.df = res.DEG.list$PCOS_D0_D2, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_DEGs_D2_18", DEG.df = res.DEG.list$PCOS_D2_D6, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_DEGs_D0_18", DEG.df = res.DEG.list$PCOS_D0_D6, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all", 
                             selected.name = "Inflammation_all", selected.genes = Inflammation.all,
                             selected.break = c(6, 15, 23), reorder.col = FALSE)

# Do heatmaps, three groups over time, Control
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[5]], x.name = "Control_all_DEGs_D0_6", DEG.df = res.DEG.list$Ctrl_D0_D2, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[5]], x.name = "Control_all_DEGs_D2_18", DEG.df = res.DEG.list$Ctrl_D2_D6, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[5]], x.name = "Control_all_DEGs_D0_18", DEG.df = res.DEG.list$Ctrl_D0_D6, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[5]], x.name = "Control_all", 
                             selected.name = "All_Selected_Genes", selected.genes = select.all.genes,
                             selected.break = c(6, 15, 23), reorder.col = FALSE) # Manually set the breaks

# Do heatmaps, three groups over time, PCOS corrected for aging DEGs
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_corr_DEGs_D0_6", DEG.df = res.DEG.list$PCOS_D0_D2_corr, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_corr_DEGs_D2_18", DEG.df = res.DEG.list$PCOS_D2_D6_corr, 
                             top.genes.n = 25, reorder.col = FALSE)
run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_corr_DEGs_D0_18", DEG.df = res.DEG.list$PCOS_D0_D6_corr, 
                             top.genes.n = 25, reorder.col = FALSE)

run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_corr", 
                             selected.name = "All_Selected_Genes", selected.genes = select.all.genes,
                             selected.break = c(6, 15, 23), reorder.col = FALSE) # Manually set the breaks

# Do heatmaps of all PCOS DEG genes
all.DEGs = res.DEG.genenames[10:12]
all.DEGs = unlist(all.DEGs)
all.DEGs = na.omit(unique(all.DEGs))
all.DEGs = unique(all.DEGs)

#run.heatmap = Do.Heatmap.fun(x.res = res.trans.list[[4]], x.name = "PCOS_all_DEGs", selected.genes = all.DEGs,
#                             selected.name = "all_DEGs",reorder.col = FALSE, convert.ensembl = TRUE)


# Plot quality control
## Generate a PCA plot
#PCA.plot = plotPCA(dds.trans, intgroup="Group") + theme_cowplot()
#ggsave2(plot = PCA.plot, filename = paste0(Output.dir, "PCA_plot_", trans.x,".pdf"))

PCA.data = plotPCA(dds.trans, intgroup="Group", returnData = TRUE)
percentVar = round(100 * attr(PCA.data, "percentVar"))
PCA.plot = ggplot(PCA.data, aes(PC1, PC2, color=Group)) +
  geom_point(size=4, stroke = 1.5, aes(shape=Group, color=Group)) +
  #scale_shape_manual(values = c(17, 15, 16, 17, 15, 16)) +
  scale_shape_manual(values = c(1, 1, 1, 16, 16, 16)) +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#E64B35FF", "#4DBBD5FF")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_cowplot()

ggsave2(plot = PCA.plot, filename = paste0(Output.dir, "PCA_plot_alt1", trans.x,".pdf"))

PCA.plot = ggplot(PCA.data, aes(PC1, PC2, color=Group)) +
  geom_point(size=4, stroke = 1.5, aes(shape=Group, color=Group)) +
  #scale_shape_manual(values = c(17, 15, 16, 17, 15, 16)) +
  scale_shape_manual(values = c(1, 1, 1, 16, 16, 16)) +
  scale_color_manual(values = rev(color.vec)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_cowplot()

ggsave2(plot = PCA.plot, filename = paste0(Output.dir, "PCA_plot_alt2", trans.x,".pdf"))

## Generate a Euclidean distance plot
sampleDists <- dist(t(assay(dds.trans)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds.trans$Group)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
euclidean.plot = pheatmap(sampleDistMatrix,
                          clustering_distance_rows=sampleDists,
                          clustering_distance_cols=sampleDists,
                          col=colors)
ggsave2(plot = euclidean.plot, filename = paste0(Supp.dir, "Euclidean_distance_plot_", trans.x,".pdf"))



