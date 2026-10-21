# Helper functions for the AmnioSC analysis

library(ggalluvial)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(viridis)
library(Polychrome)
library(cowplot)
library(patchwork)
library(egg)
library(SingleCellExperiment)
library(scran)
library(Seurat)
library(SeuratObject)

# PAPER SAMPLE AND DONOR NAMES -------------------------------------------------

SAMPLE.PAPER <- c("AF_1a", "AF_1b", "AF_2", "AF_3", "AF_4", "AF_5", "AF_6", "AF_7")
names(SAMPLE.PAPER) <- c("1AFa", "1AFb", "2AF", "4AF", "5AF", "6AF", "7AF", "8AF")

DONOR.PAPER <- c("AF_1", "AF_2", "AF_3", "AF_4", "AF_5", "AF_6", "AF_7")
names(DONOR.PAPER) <- c("1AF", "2AF", "4AF", "5AF", "6AF", "7AF", "8AF")

# HELPER FUNCTIONS -------------------------------------------------------------

calcMDS <- function(pseudobulk, max_dim, ID) {
  df <- data.frame(ID = rownames(pseudobulk$samples))
  rownames(df) <- df$ID
  colnames(df) <- ID
  for (i in 1:(max_dim - 1)) {
    currMDS <- limma::plotMDS(edgeR::cpm(pseudobulk, log=TRUE), 
                              dim.plot = c(i, i + 1), 
                              gene.selection = "common")
    df[, paste0("Dim_", i)] <- currMDS$x
    if (i == max_dim - 1) {
      df[, paste0("Dim_", max_dim)] <- currMDS$y
    }
  }
  return(df)
}

plotMDS.Dims1_5 <- function(MDS.object, color_var, COLOURS = NULL, makeFactor = FALSE) {
  if (makeFactor) {
    MDS.object[, color_var] <- factor(MDS.object[, color_var])
  }
  
  p1 <- ggplot(MDS.object, aes(x = Dim_1, y = Dim_2, color = get(color_var))) +
    geom_point(size = 5) +
    theme_cowplot(12) +
    no_legend
  
  p2 <- ggplot(MDS.object, aes(x = Dim_2, y = Dim_3, color = get(color_var))) +
    geom_point(size = 5) +
    theme_cowplot(12) +
    no_legend
  
  p3 <- ggplot(MDS.object, aes(x = Dim_3, y = Dim_4, color = get(color_var))) +
    geom_point(size = 5) +
    theme_cowplot(12) +
    bottom_legend +
    labs(color = color_var)
  
  p4 <- ggplot(MDS.object, aes(x = Dim_4, y = Dim_5, color = get(color_var))) +
    geom_point(size = 5) +
    theme_cowplot(12) +
    no_legend
  
  if (!is.null(COLOURS)) {
    p1 <- p1 + scale_color_manual(values = COLOURS[levels(MDS.object[, color_var])])
    p2 <- p2 + scale_color_manual(values = COLOURS[levels(MDS.object[, color_var])])
    p3 <- p3 + scale_color_manual(values = COLOURS[levels(MDS.object[, color_var])])
    p4 <- p4 + scale_color_manual(values = COLOURS[levels(MDS.object[, color_var])])
  }
  
  if (is.null(COLOURS) & !makeFactor) {
    p1 <- p1 + scale_color_viridis(trans = "log10")
    p2 <- p2 + scale_color_viridis(trans = "log10")
    p3 <- p3 + scale_color_viridis(trans = "log10")
    p4 <- p4 + scale_color_viridis(trans = "log10")
  }
  
  print((p1 + p2) / (p3 + p4))
}

createAlluvialLongDF <- function(names, label1, label2, label3=NULL, label4=NULL, label5=NULL) {
  # Creates and returns a long format data.frame for alluvial plots
  # Accepts a minimum of 3 vectors (name and two labels); 3 additional labels are accepted
  alluvial_change_df <- data.frame(names = names,
                                   label1 = label1,
                                   label2 = label2)
  
  alluvial_change_df$label3 <- label3
  alluvial_change_df$label4 <- label4
  alluvial_change_df$label5 <- label5
  
  rownames(alluvial_change_df) <- alluvial_change_df$names
  
  alluvial_change_df_long <- melt(alluvial_change_df, id.vars = "names")
  colnames(alluvial_change_df_long) <- c("names", "label", "label_type")
  
  alluvial_change_df_long$label_type <- factor(alluvial_change_df_long$label_type,
                                               levels = names(sort(
                                                 table(alluvial_change_df_long$label_type), 
                                                 decreasing = TRUE)))
  return(alluvial_change_df_long)
}

lognormCounts <- function(x, log = TRUE, prior.count = 0.5)
  # Function to log normalise a counts matrix to median library size
  # Input is counts matrix
  # Genes as rows, cells as columns
  # Belinda Phipson
  # 26 February 2020
  # https://github.com/bphipson/Human_Development_snRNAseq/blob/master/code/normCounts.R
{
  x <- as.matrix(x)
  lib.size <- colSums(x)
  M <- median(lib.size)
  if(log){
    prior.count.scaled <- lib.size/mean(lib.size)*prior.count
    lib.size <- lib.size + 2*prior.count.scaled
    log2(t((t(x)+prior.count.scaled)/lib.size*M))
  }
  else t(t(x)/lib.size*M)
}

# order factor by most common labels
sortedFactor <- function(values, decreasing = TRUE) {
  # Accepts a vector of values/column from dataframe, etc.
  # Converts to factor with levels sorted by most common value
  sorted_values <- factor(values,
                          levels = names(sort(table(values), 
                                              decreasing = decreasing)))
  return(sorted_values)
}

integrateDataset <- function(seu, split_var = "Dataset", do_split = TRUE, 
                             k.weight=100, npcs=50) {
  # Accepts a Seurat object containing multiple datasets and performs data integration
  # Can alternatively accept a list of Seurat objects if `do_split` is TRUE
  # Based upon the integration approach by Jovana Maksimovic
  
  if (do_split == TRUE) {
    DefaultAssay(seu) <- "RNA"
    seu_list <- SplitObject(seu, split.by = split_var)
  } else {
    seu_list <- seu
  }
  
  seu_list <- lapply(X = seu_list, FUN = SCTransform, method = "glmGamPoi")
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = seu_list,
                                        nfeatures = 3000)
  seu_list <- PrepSCTIntegration(object.list = seu_list,
                                 anchor.features = features)
  seu_list <- lapply(X = seu_list, FUN = RunPCA, features = features, npcs = npcs)
  
  anchors <- FindIntegrationAnchors(object.list = seu_list, 
                                    anchor.features = features,
                                    normalization.method = "SCT",
                                    reduction = "rpca",
                                    dims = 1:30,
                                    k.anchor = 20)
  # creates an 'integrated' data assay
  seu_int <- IntegrateData(anchorset = anchors,
                           k.weight = k.weight,
                           normalization.method = "SCT",
                           dims = 1:30)
  DefaultAssay(seu_int) <- "integrated"
  
  seu_int <- RunPCA(seu_int, npcs = 30, verbose = FALSE)
  seu_int <- RunUMAP(seu_int, dims = 1:30, verbose = FALSE)
  
  return(seu_int)
}

findMarkerGenes <- function(counts, group, samples=NULL) {
  # Accepts a counts matrix and associated group and sample (optional) vectors
  # Calculates logcounts using lognormCounts() function
  # Performs marker gene analysis by group using limma
  # Returns the fit.cont object from eBayes
  
  y.org <- DGEList(counts)
  logcounts <- lognormCounts(y.org, log = TRUE, prior.count = 0.5)
  
  # create design matrix
  if (is.null(samples)) {
    design <- model.matrix(~ 0 + group)
  } else {
    design <- model.matrix(~ 0 + group + samples)
  }
  colnames(design)[1:(length(levels(group)))] <- levels(group)
  
  # create contrast matrix
  mycont <- matrix(NA, ncol = length(levels(group)), 
                   nrow = length(levels(group)))
  rownames(mycont) <- colnames(mycont) <- levels(group)
  diag(mycont) <- 1
  mycont[upper.tri(mycont)] <- -1/(length(levels(factor(group))) - 1)
  mycont[lower.tri(mycont)] <- -1/(length(levels(factor(group))) - 1)
  
  # fill out remaining rows with 0s
  zero.rows <- matrix(0, ncol = length(levels(group)),
                      nrow = (ncol(design) - length(levels(group))))
  fullcont <- rbind(mycont, zero.rows)
  rownames(fullcont) <- colnames(design)
  
  fit <- lmFit(logcounts, design)
  fit.cont <- contrasts.fit(fit, contrasts = fullcont)
  fit.cont <- eBayes(fit.cont, trend = TRUE, robust = TRUE)
  
  return(fit.cont)
}

extractMarkerGenes <- function (treat.all, group_names, n_markers=10,
                                group_name_prefix="") {
  # Accepts treat.all object from treat() and associated group names
  # Extracts the desired number of markers from the object per group and 
  # returns these as a named list
  # Group name prefix allows for numerical group names to be altered to allow
  # them to be used as dataframe column names, etc. easily
  
  markers_topN <- list()
  for (i in c(1:length(group_names))) {
    top <- topTreat(treat.all, coef = i, n = Inf, sort.by = "t")
    markers_topN[[i]] <- rownames(top)[top$logFC>0][1:n_markers]
  }
  names(markers_topN) <- paste0(group_name_prefix, group_names)
  return(markers_topN)
}

performGO <- function(markers, group_names, group_name_prefix="",
                                 n_go=10) {
  # Accepts a named list of markers for each group and associated group names
  # Returns a list containing the top N GO terms for each group
  # Group name prefix allows for numerical group names to be altered to allow
  # them to be used as dataframe column names, etc. easily
  
  go_terms <- list()
  for (i in c(1:length(group_names))) {
    curr_markers <- markers[[i]]
    
    # convert symbol to entrez gene names
    curr_markers <- mapIds(org.Hs.eg.db, keys = curr_markers,
                           column = "ENTREZID", keytype = "SYMBOL")
    
    # Gene Ontology
    curr_markers_go <- goana(unique(curr_markers), species = "Hs", universe = NULL)
    curr_markers_go <- curr_markers_go[order(curr_markers_go$P.DE),]
    # Only keeping biological process terms that are not overly general.
    curr_markers_go <- curr_markers_go[curr_markers_go$Ont=="BP",]
    curr_markers_go <- curr_markers_go[curr_markers_go$N <= 200,]
    go_terms[[paste0(group_name_prefix, group_names[i])]] <- curr_markers_go[1:n_go, 1]
    
  }
  return(go_terms)
}

createCelltypistResultsDF <- function(predicted_labels, probability_matrix, decision_matrix, ref_name="") {
  # Accepts output from Celltypist annotation tool and reformats output so it 
  # can be added as Seurat/SCE metadata 
  
  # reformat predicted labels data frame
  rownames(predicted_labels) <- predicted_labels$X
  
  # replace problematic characters with full stops
  predicted_labels$predicted_labels <- gsub("[^[:alnum:]_]", ".", predicted_labels$predicted_labels)
  predicted_labels$majority_voting <- gsub("[^[:alnum:]_]", ".", predicted_labels$majority_voting)
  
  predicted_labels$prob_pred <- 0
  predicted_labels$decision_pred <- 0
  predicted_labels$prob_vote <- 0
  predicted_labels$decision_vote <- 0
  
  # add rownames to probability and decision matrices
  rownames(probability_matrix) <- probability_matrix$X
  rownames(decision_matrix) <- decision_matrix$X
  
  # add in the corresponding probability and decision value for each of the assigned labels
  for (cell in rownames(predicted_labels)) {
    predicted_labels[cell, "prob_pred"] <- probability_matrix[cell, predicted_labels[cell, "predicted_labels"]]
    predicted_labels[cell, "prob_vote"] <- probability_matrix[cell, predicted_labels[cell, "majority_voting"]]
    
    predicted_labels[cell, "decision_pred"] <- decision_matrix[cell, predicted_labels[cell, "predicted_labels"]]
    predicted_labels[cell, "decision_vote"] <- decision_matrix[cell, predicted_labels[cell, "majority_voting"]]
  }
  
  # remove Roman numerals from end of cell type names
  predicted_labels$predicted_labels <- sub("_[VXI]*$", "", predicted_labels$predicted_labels)
  predicted_labels$majority_voting <- sub("_[VXI]*$", "", predicted_labels$majority_voting)
  predicted_labels$predicted_labels <- gsub(" [0-9]", "", predicted_labels$predicted_labels)
  predicted_labels$majority_voting <- gsub(" [0-9]", "", predicted_labels$majority_voting)
  
  # remove X added to start of cell names if added
  predicted_labels$X <- sub("^X", "", predicted_labels$X)
  # replace . in cell name with -
  predicted_labels$X <- gsub("\\.", "-", predicted_labels$X)
  # update rownames with altered cell names
  rownames(predicted_labels) <- predicted_labels$X
  
  predicted_labels <- predicted_labels[, c("predicted_labels", "majority_voting",
                                           "prob_pred", "prob_vote", 
                                           "decision_pred", "decision_vote")]
  
  colnames(predicted_labels) <- paste0(ref_name, ".", colnames(predicted_labels))
  
  return(predicted_labels)
}

plotQCViolin <- function(seu_object, x, colour, COLOURS, title = NULL) {
  seu_object@meta.data %>%
    ggplot(aes(x = get(x),
               y = nCount_RNA)) +
    geom_violin(scale = "area") +
    geom_quasirandom(aes(color = get(colour))) +
    scale_y_log10() +
    theme_cowplot(12) +
    xlab(x) +
    scale_colour_manual(values = COLOURS[levels(seu_object@meta.data[, colour])]) +
    NoLegend() -> p1
  
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title)
  }
  
  seu_object@meta.data %>%
    ggplot(aes(x = get(x),
               y = nFeature_RNA)) +
    geom_violin(scale = "area") +
    geom_quasirandom(aes(color = get(colour))) +
    scale_y_log10() +
    theme_cowplot(12) +
    xlab(x) +
    labs(color = colour) + 
    scale_colour_manual(values = COLOURS[levels(seu_object@meta.data[, colour])]) -> p2
  
  print(p1 | p2)
}

# Drop any markers not detected in this dataset and flag them for visibility,
# rather than letting DotPlot() error out on a missing feature.
filterPresentGenes <- function(genes, object, assay = "RNA") {
  available <- genes[genes %in% rownames(object[[assay]])]
  missing <- setdiff(genes, available)
  if (length(missing) > 0) {
    message("Not detected in dataset, skipping: ", paste(missing, collapse = ", "))
  }
  return(available)
}

generateGroupProportionsDF <- function(seu_int, group_var, donor_var = "Donor") {
  # Prepare a counts dataframe for proportion plots
  counts_df <- seu_int@meta.data %>%
    group_by(.data[[donor_var]], .data[[group_var]]) %>%
    summarise(
      n = n(),
      early_late_gestation = unique(early_late_gestation),
      .groups = "drop"
    ) %>%
    # Fill in zero counts for any donor x cluster combinations that are absent
    complete(.data[[donor_var]], .data[[group_var]], fill = list(n = 0)) %>%
    # Re-fill early_late_gestation for the completed rows
    group_by(.data[[donor_var]]) %>%
    fill(early_late_gestation, .direction = "downup") %>%
    ungroup() %>%
    mutate(
      # Set Early as the reference level so the model estimates Late vs Early
      early_late_gestation = factor(early_late_gestation, levels = c("Early", "Late"))
    )
  
  counts_df <- counts_df %>% group_by(.data[[donor_var]]) %>% mutate(total_Donor = sum(n))
  counts_df$proportion_cluster <- counts_df$n / counts_df$total_Donor
  return(counts_df)
}

# PRETTY PLOTS FOR THE MANUSCRIPT ----------------------------------------------

plotUMAP.Paper <- function(seu_int, color_by, colour_scheme, title = "", subtitle="",
                           label_size = 8, label_clusters = TRUE) {
  # UMAP plot formatted for the paper
  umap <- DimPlot(seu_int, reduction = 'umap', label = FALSE, 
                  group.by = color_by, pt.size = 2) + theme_void(12) +
    scale_color_manual(values = colour_scheme) + blank_axes +
    no_legend +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    theme(plot.title = element_markdown(face = "bold", size = 24, 
                                        hjust = 0.5, family = "arial"),
          plot.subtitle = element_markdown(face = "italic", size = 22, 
                                           hjust = 0.5, family = "arial"))
  
  if (isTRUE(label_clusters)) {
  umap <- LabelClusters(umap,
                        id = color_by, fontface = "bold",
                        color = "black", size = label_size, max.overlaps = 100,
                        repel = TRUE, box = FALSE, fill = NA,
                        alpha = 1, seed = 1234, family = "arial")
  }
  
  return(umap)
}

plotNCells.Paper <- function(seu_int, x_var, color_by, colour_scheme, title = "", subtitle="",
                             xlab = "", ylab = "", ymax = NULL) {
  seu_int@meta.data$x_plot <- sortedFactor(seu_int@meta.data[,x_var],
                                           decreasing = FALSE)
  seu_int@meta.data %>%
    ggplot(aes(x = x_plot, fill = get(color_by))) +
    geom_bar() +
    geom_text(stat='count', aes(label=..count..), 
              hjust = -0.2, size = 6, family = "arial") +
    theme_linedraw(12) +
    scale_fill_manual(values = colour_scheme) +
    no_legend + coord_flip() +
    ylab(ylab) + xlab(xlab) + labs(title = title, subtitle = subtitle) +
    theme(plot.title = element_markdown(face = "bold", size = 24, family = "arial"),
          plot.subtitle = element_markdown(face = "italic", size = 22, family = "arial"),
          axis.title = element_text(size = 22, family = "arial", face = "bold"),
          axis.text = element_text(size = 18, family = "arial")) -> bar_plot
  
  if (!is.null(ymax)) {
    bar_plot <- bar_plot + ylim(0, ymax)
  }
  
  return(bar_plot)
}

plotProportion.Paper <- function(plot_df, curr_cluster, title = "", 
                                 sccomp_result, font_scaling = 0) {
  sccomp_result <- sccomp_result %>% filter(cluster_names == curr_cluster,
                                            parameter == "early_late_gestationLate")
  
  ggplot(plot_df %>% filter(cluster_names == curr_cluster), 
         aes(x = early_late_gestation, y = proportion_cluster, 
             color = early_late_gestation)) +
    geom_jitter(width = 0.15, size = 4) +
    theme_linedraw(16) + no_legend + 
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_color_manual(values = ALL_GEST_COLOURS[-2]) +
    ylab("") + xlab("") +
    labs(title = title,
         subtitle = paste0("FDR: ", 
                           round(sccomp_result$c_FDR, 6))) +
    theme(plot.title = element_markdown(face = "bold", size = 18 + font_scaling, family = "arial"), 
          plot.subtitle = element_text(face = "italic", size = 16 + font_scaling, family = "arial"),
          axis.title = element_text(size = 16 + font_scaling, family = "arial", face = "bold"),
          axis.text = element_text(size = 14 + font_scaling, family = "arial")) -> proportion_plot
  
  return(proportion_plot)
}

# PLOTTING SHORTCUTS -----------------------------------------------------------

# shortcut to rotate x-axis labels on ggplot2 plots
rotate_x_labels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                    hjust=1))

rotate_x_labels_45 <- theme(axis.text.x = element_text(angle = 45, hjust=1))

# shortcut to remove y-axis labels on ggplot2 plots
no_y_labels <- theme(axis.title.y = element_blank(), 
                     axis.text.y = element_blank())

# shortcut to remove legend on ggplot2 plots
no_legend <- theme(legend.position = "none")

# shortcut to place legend at bottom of ggplot2 plots
bottom_legend <- theme(legend.position = "bottom")

# remove the axis labels and ticks for both axes
blank_axes <- theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
