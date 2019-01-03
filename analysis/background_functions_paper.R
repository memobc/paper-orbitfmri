###############################################################################################################

# R functions to summarise and plot background connectivity data from CONN:

# Rose Cooper - last updated Jan 2019

###############################################################################################################


## 1. Format conn data ------------------------------------------------------------------------------------ #
format_conn <- function(myDir, seeds, subjects, model, event){
  
  # loop through subjects to extract ROI-to-ROI matrix fron CONN first level analyses
  connMatrix = array(data = NA, c(length(seeds),length(seeds),length(subjects)))
  for (s in 1:length(subjects)) {

    # grab fisher z/beta connectivity matrix from CONN output:
    subMatrix <- data.matrix(read.csv(paste(myDir,model,'/',subjects[s],'/ROI-ROI_',event,'.csv',sep = ""), header = TRUE))
    # set matrix diagonal to zero
    diag(subMatrix) <- 0
    connMatrix[,,s] <- subMatrix
    
  } #end of loop through subjects
  
  return(connMatrix)
}
# ------------------------------------------------------------------------------------------------ #



## 2. Get corrected p-values per ROI Pair -------------------------------------------- #
get_conn_pValues <- function(connMatrix, seeds, sOrder, type, th){
  
  ### INPUT:
  # connMatrix  --> a matrix of size seed x target x subject containing the z/r estimates
  # seeds       --> a vector of seed names used for the conn analysis.
  # sOrder      --> order to sort ROIs
  # type        --> type of t-test - 'two-tailed', 'greater', 'less'
  # th          --> r threshold to define 'connections'
  
  ### OUTPUT:
  # connPvalues  --> a data.frame with the FDR corrected p value per seed-target comparison and * to indicate < .05
  
  connPVals = array(data = 1, c(nrow(connMatrix),ncol(connMatrix)))
  for (seed in 1:nrow(connMatrix)) {
    for (targ in 1:ncol(connMatrix)) {
      if (seed != targ) {
        connVec <- connMatrix[seed,targ,] #seed-target for all subjects
        if (type == 'greater') {
          # one-sample t-test for this ROI-ROI connection across subjects
          connPVals[seed,targ] <- t.test(connVec, mu=fisherz(th), alternative = type)$p.value
        } else if (type == 'two.sided') {  #difference between conditions
          connPVals[seed,targ] <- t.test(connVec, mu=0, alternative = type)$p.value
        }
      }
    } # end of loop through targets
  }  # end of loop through seeds
  
  # select half -- symmetrical
  connPVals[lower.tri(connPVals, diag = FALSE)]<- NA
  
  ### b) adjust p values for multiple comparisons based on all possible connections
  ### FDR adjustment:
  connPVals <- data.frame(connPVals)
  colnames(connPVals) <- seeds
  connPVals$seeds <-seeds
  connPVals$seeds <- as.factor(connPVals$seeds)
  connPVals$seeds = factor(connPVals$seeds,levels(connPVals$seeds)[sOrder])
  melted_connPVals <- melt(connPVals, id="seeds", na.rm = TRUE)
  
  # adjust
  melted_connPVals[melted_connPVals == 1] <- NA #remove diagonal so not included in multiple comparison correction
  melted_connPAdjust <- p.adjust(melted_connPVals$value, method = "fdr", n = (length(melted_connPVals$value) - length(seeds)))
  melted_connPVals$value <- melted_connPAdjust #replace original p values with adjusted
  
  # add significance asterix for overlay on connection matrix:
  melted_connPVals$sig <- ''
  for (r in 1:nrow(melted_connPVals)) {
    if (!is.na(melted_connPVals$value[r])) { 
      if (as.numeric(melted_connPVals$value[r]) < 0.05) {  #if significant FDR corrected, add asterix
        melted_connPVals$sig[r] <- '*'
      } 
    }
  }  # end of loop through seeds 
  
  melted_connPVals <- data.frame(melted_connPVals)
  return(melted_connPVals)
}
# ------------------------------------------------------------------------------------------------ #



## 3. Plot mean (group-level) conn matrix --------------------------------------------------------- #
plot_meanconn <- function(connMatrix, melted_connPVals, seeds, sOrder, title, type, th){
  
  ### INPUT:
  # connMatrix  --> a matrix of size seed x target x subject containing the z/r estimates
  # melted_connPVals --> data frame of p values and * from get_conn_pValues
  # seeds      --> a vector of seed names that were those used for the conn analysis.
  # sOrder    --> order to sort ROI levels
  # title     --> title for plot
  # type      --> type of t-test - 'two-tailed', 'greater', 'less'
  # th        --> r threshold to define 'connections'
  
  ### OUTPUT:
  # geom_tile ggplot
  
  ##### Calculate mean matrix across subjects (already fisher z scores)
  meanconn <- apply(connMatrix, c(1,2), mean)

  #convert to r to visualize
  meanconn <- fisherz2r(meanconn)

  ## symmetrical - select only half
  meanconn[lower.tri(meanconn, diag = FALSE)]<- NA
  
  meanconn <- data.frame(meanconn)
  colnames(meanconn) <- seeds
  meanconn$seeds <-seeds
  meanconn$seeds <- as.factor(meanconn$seeds)
  meanconn$seeds = factor(meanconn$seeds,levels(meanconn$seeds)[sOrder])
  melted_meanconn <- melt(meanconn, id="seeds", na.rm = TRUE)
  melted_meanconn$value[melted_meanconn$value == 0] <- NA
  
  # now r values are same format as p values, set non-sig correlations to NA
  melted_meanconn$value[melted_connPVals$value >= 0.05] <-NA

  # plot  
  myCol <- c("dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2",
             "mediumorchid","mediumorchid",
             "firebrick2","firebrick2","firebrick2","firebrick2","firebrick2")
  
  # set scale for plot
  if (type == 'greater') {
     limits = c(th,0.8)
  } else if (type == 'two.sided') {  #difference between conditions
    limits = c(-0.25,0.25)
  }

  p <- ggplot(data = melted_meanconn, aes(x=seeds, y=variable, fill=value)) + 
    geom_tile(color = "white")+
    ggtitle(title) +
    scale_fill_gradient2(low = "#99d594", mid = "#fee08b", high = "#d53e4f", na.value = "grey96",
                         midpoint = median(limits), limit = limits, name="Connectivity") +
    xlab('ROI') + ylab('ROI') +
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=24, family="Helvetica", colour=myCol),
          axis.text.y = element_text(angle = 45, hjust=1, size=24, family="Helvetica", colour=myCol),
          axis.title = element_text(size=28,family="Helvetica"), 
          legend.direction = "horizontal", plot.title = element_text(hjust = 0.5, size=32),
          legend.title = element_text(size = 15), legend.key.size = unit(1,"cm"),
          legend.text = element_text(size = 12), legend.position = "bottom",
          panel.background = element_blank()) + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
  
  return(p)
}
# ------------------------------------------------------------------------------------------------ #



# 4. Run modularity on connectivity matrix ------------------------------------------------------- #
run_modularity <- function(connMatrixE, connMatrixR, subjects, th){
  
  # Runs Network Toolbox louvain function on a correlation matrix (with diagnonal set to 0)
  # returns the modularity value (Q) per task
  
  Q <- data.frame(array(data = NA, c(length(subjects)*2,3)))
  colnames(Q) = c('Subject','Task','Modularity')
  row = 0
  for (sub in 1:length(subjects)) {
    for (t in 1:2) {
      if (t == 1){
        curMatrix <- fisherz2r(connMatrixE[,,sub])  #convert back to r - CONN saves as z for later averaging
        task = 'Encoding'
      } else {
        curMatrix <- fisherz2r(connMatrixR[,,sub])  #convert back to r - CONN saves as z for later averaging
        task = 'Retrieval'
      }
    row = row+1
    # remove all low (below-threshold) correlations
    curMatrix[curMatrix < th] <- 0
    # compute modularity
    modularity <- louvain(curMatrix)
    # add info
    Q$Subject[row] <- sub
    Q$Task[row] <- task
    Q$Modularity[row] <- modularity$Q
    }
  }
  
  return(Q)
  
}
# ------------------------------------------------------------------------------------------------ #



# 5. Run within vs between strength on connectivity matrix ------------------------------------------------------- #
run_network<- function(connMatrixE, connMatrixR, subjects, networks, th){
  
  Q <- data.frame(array(data = NA, c(length(subjects)*4,4)))
  colnames(Q) = c('Subject','Task','Measure','Strength')
  row = 0
  for (sub in 1:length(subjects)) {
    for (t in 1:2) {
      if (t == 1){
        curMatrix <- fisherz2r(connMatrixE[,,sub])  #convert back to r - CONN saves as z for later averaging, density calculated on correlation coefficients
        task = 'Encoding'
      } else {
        curMatrix <- fisherz2r(connMatrixR[,,sub])  #convert back to r - CONN saves as z for later averaging
        task = 'Retrieval'
      }

      # remove all low (below-threshold) correlations before compute density
      curMatrix[curMatrix < th] <- 0
      diag(curMatrix) <- NA

      # get all within and between network connections:
      within <- NA
      between <- NA
      for (n in 1:length(unique(networks))) {
        within <- c(within, curMatrix[networks == n, networks == n])
        between <- c(between, curMatrix[networks == n, networks != n])
      }

      row = row+1
      Q$Subject[row] <- sub
      Q$Task[row] <- task
      Q$Measure[row] <- 'Within'
      Q$Strength[row] <- mean(within, na.rm = TRUE)
      
      row = row+1
      Q$Subject[row] <- sub
      Q$Task[row] <- task
      Q$Measure[row] <- 'Between'
      Q$Strength[row] <- mean(between, na.rm = TRUE)
    }
  }
  
  return(Q)
}
# ------------------------------------------------------------------------------------------------ #