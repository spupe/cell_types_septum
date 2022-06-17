
#### "Septum_Functions"

##### Medial Septum sequencing project 
##### Script by Stefano Pupe
  
##### VISUALIZATION - modify_vlnplot, StackedVlnPlot
# VlnPlot from Seurat kind of sucks for making figures, so I'm using StackedVlnPlot(), written by Tommy Tang (slightly modified)

# This function modifies Seurat's VlnPlot to make it prettier
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0,
                          title = "",
                          plot.margin = margin(0,0,0,0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle(title) + 
    theme(legend.position = "none", 
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1),angle=0),  
          axis.text.y = element_blank(), 
          plot.margin = plot.margin ) 
  return(p)
}


# This small function is required as StackedVlnPlot needs to know the y limit for all the panels 
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

# Main function, it creates a patchwork composed of violin plots (produced by the modify_vlnplot function) and stacks them nicely
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          verbose = FALSE,
                          
                          ...) {
  
  plot_list<- purrr::map(rev(features), function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle=0), axis.ticks.x = element_line()) +annotate(geom = 'text', 
                                                                                      fontface = 'italic')
  
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x + scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  
  return(p)
}

# VlnPlot does not plot genes that are not found in the list. The problem is worse with StackedVlnPlot, which doesn't even run
# if there are any missing genes. We have a use case where we will pass a large number of genes to StackedVlnPlot (Figure S3),
# so we need this quick function

checkGenes <- function(seuratObj, 
                       markerstocompare,
                       verbose = FALSE,
                       ...) {
  gene_list <- rownames(seuratObj@assays$RNA)
  genesfound <- intersect(gene_list,markerstocompare)
  genesnotfound <- setdiff(markerstocompare,gene_list)
  
  output <- list(genesfound,genesnotfound)
  
  if (length(output[[2]]>0)) {
     print(paste0("Gene not found in Seurat object: ", output[[2]]))
  }
  
  return(output)
} 


#### CLUSTER TESTING - SetIfNull, RFClassifier, AssessSplit, ClassifyCells
# SetIfNull, RFClassifier and AssessSplit functions from Seurat v2, adapted to v3 by me. 
# The latter two functions are useful for evaluating the robustness of clustering, and ClassifyCells is sort of a frontend for them

# Very simple function, simply checks if x is null and if so, replaces it with a default value
SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}

# Function extracted from Seurat v2, very slighly changed due to "grammar" updates in seurat object handling in v3
# The function takes the Random Forest implementation from the 'ranger' package, and uses it on a seurat object, trying to see
# if it can accurately predict the identity (class) based on the expression of a set of genes (training.genes)

 
BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
  
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(
    x = training.genes,
    default = rownames(x = object@assays$integrated@data)
  )
  training.data <- as.data.frame(
    x = t(
      x = as.matrix(
        x = object@assays$integrated@data[training.genes, ]
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    message("Training Classifier ...")
  }
  classifier <- ranger::ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

# This function is a specific application of the one before, also from Seurat v2. It compares two specific idents by using 
# BuildRFClassifier, and reporting the prediction error (out of bag error)

AssessSplit <- function(
  object,
  node,
  cluster1,
  cluster2,
  print.output = TRUE,
  verbose = TRUE,
  ...
) {
  tree <- object@tools$BuildClusterTree[[1]]
  if (! missing(x = node)){
    if (! missing(x = cluster1) || ! missing(x = cluster2)) {
      warning("Both node and cluster IDs provided. Defaulting to using node ID")
    }
    possible.nodes <- c(
      DFT(tree = tree, node = tree$edge[,1][1]),
      tree$edge[,1][1]
    )
    if (! node %in% possible.nodes) {
      stop("Not a valid node")
    }
    split <- tree$edge[which(x = tree$edge[,1] == node), ][,2]
    group1 <- DFT(tree = tree, node = split[1], only.children = TRUE)
    group2 <- DFT(tree = tree, node = split[2], only.children = TRUE)
    if (any(is.na(x = group1))) {
      group1 <- split[1]
    }
    if (any(is.na(x = group2))) {
      group2 <- split[2]
    }
  } else {
    group1 <- cluster1
    group2 <- cluster2
  }
  group1.cells <- WhichCells(object = object, ident = group1)
  group2.cells <- WhichCells(object = object, ident = group2)
  assess.data <- subset(
    x = object,
    cells = c(group1.cells, group2.cells)
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group1.cells,
    value = "g1"
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group2.cells,
    value = "g2"
  )
  rfc <- BuildRFClassifier(
    object = assess.data,
    training.genes = assess.data@assays$integrated@var.features,
    training.classes = assess.data@active.ident,
    ...
  )
  oobe <- rfc$prediction.error
  if (print.output) {
    print(paste0("Out of Bag Error: ", round(x = oobe, digits = 4) * 100, "%"))
  }
  #return(oobe) 
  return(rfc)
}

# Another Seurat v2 function, it is a frontend to BuildRFclassifier. It takes in a gene set and a list of 
# idents, and gives out a specific prediction about which ident it thinks the cells expressing those genes
# belong to

ClassifyCells <- function(
  object,
  classifier,
  training.genes = NULL,
  training.classes = NULL,
  new.data = NULL,
  ...
) {
  
  if (missing(classifier)){
    classifier <- BuildRFClassifier(
      object = object,
      training.genes = training.genes,
      training.classes = training.classes,
      ...
    )
  }
  
  features <- classifier$forest$independent.variable.names
  genes.to.add <- setdiff(x = features, y = rownames(x = new.data))
  data.to.add <- matrix(
    data = 0,
    nrow = length(x = genes.to.add),
    ncol = ncol(x = new.data)
  )
  rownames(x = data.to.add) <- genes.to.add
  new.data <- rbind(new.data, data.to.add)
  new.data <- new.data[features, ]
  new.data <- as.matrix(x = t(x = new.data))
  message("Running Classifier ...")
  prediction <- predict(classifier, new.data)
  new.classes <- prediction$predictions
  return(new.classes)
}

# rfClassify function built by Mark Cembrowski, adapted by me
# This is in turn a frontend to classifyCells, it repeatedly (for the number of its) compares the predictions
# generated by the random forest classifier, with a growing proportion of the dataset, and then plots the accuracy of 
# the classifier according to the number of cells that was used to train it


rfClassify <- function(seuratObj,its=100){ 
  nCells <- length(colnames(seuratObj))
  nTrains <- c(round((nCells/100)*5),round((nCells/100)*10),round((nCells/100)*25),round((nCells/100)*50),round((nCells/100)*75)) # This is the percentage of cells in the training dataset, change this if interested in more granularity
  perf <- matrix(nrow=length(nTrains),ncol=its)
  nCells <- seuratObj@assays$integrated@data@Dim[2]
  for (jj in 1:length(nTrains)){
    nTrain <- nTrains[jj]
    for (ii in 1:its){
      set.seed((jj-1)*its+ii)
      indTrain <- sample(1:nCells,nTrain)
      indTest <- 1:nCells
      indTest <- indTest[!indTest%in%indTrain]
      
      train <- subset(seuratObj, cells = seuratObj@assays$integrated@data@Dimnames[[2]][indTrain])
      test <- subset(seuratObj, cells = seuratObj@assays$integrated@data@Dimnames[[2]][indTest])
      classesPred <- ClassifyCells(
        object = train,
        training.classes = train@active.ident,
        new.data = as.matrix(test@assays$integrated@data)
      )
      diff <- as.character(classesPred)!=as.character(test@active.ident)
      
      successRat <- round(sum(abs(diff)<0.1)/length(diff)*100)
      
      perf[jj,ii] <- successRat
      print(paste('Correct class predicted with following percent:',successRat))
    }
  }
  rownames(perf) <- nTrains
  
  # Do stats
  df <- data.frame(
    mu=apply(perf,1,mean),
    sd=apply(perf,1,sd),
    lo=apply(perf,1,mean)-apply(perf,1,sd),
    hi=apply(perf,1,mean)+apply(perf,1,sd),
    nTrain=nTrains
  )
  
  # Plot
  gg <- ggplot(df,aes(x=nTrain,y=mu))
  gg <- gg + geom_ribbon(aes(ymin=lo,ymax=hi),colour='grey',fill='grey')
  gg <- gg + geom_line() + geom_point() + theme_bw()
  gg <- gg + xlab('Number in training dataset')
  gg <- gg + ylab('Percent success')
  gg <- gg + coord_cartesian(xlim=c(0,nCells),ylim=c(0,100))	
  print(gg)	
  
  invisible(df)
}

#classify_result <- rfClassify(septum.integrated) # Run this for testing the function; note that it could take quite long (up to 1h) for full run
