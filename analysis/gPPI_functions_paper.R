###############################################################################################################

# Functions to format and test gPPI first level data in R

# Rose Cooper - last updated Jan 2019

###############################################################################################################

## 1. Format conn data ------------------------------------------------------------------------------------ #
format_conn <- function(myDir, seeds, subjects, model, contrast, event){
  
  # loop through subjects to extract ROI-to-ROI matrix fron CONN first level analyses
  connMatrix = array(data = NA, c(length(seeds),length(seeds),length(subjects)))
  for (s in 1:length(subjects)) {
    
    # grab fisher z/beta connectivity matrix from CONN output:
    subMatrix <- data.matrix(read.csv(paste(myDir,model,'_ROI-to-ROI/',subjects[s],'/ROI-ROI_',
                                            event,'x',contrast,'^1','.csv',sep = ""), header = TRUE))
    # set matrix diagonal to zero
    diag(subMatrix) <- 0
    connMatrix[,,s] <- subMatrix
    
  } #end of loop through subjects
  
  return(connMatrix)
}
# ------------------------------------------------------------------------------------------------ #



## 2. Overall change in connectivity for aHipp and pHipp (seed/target) --------------------------- #
mean_hipp_change<- function(subjects, ppiMatrix, seeds, networks){
  
  meanConn <- data.frame(array(NA, c(length(subjects)*4,4)))
  colnames(meanConn) <- c('Subject','Network','Region','MeanConnectivity')
  row=0
  for (s in 1:length(subjects)) {
    curMatrix <- ppiMatrix[,,s]
    
    for (hipp in c('aHIPP','pHIPP')) {
      for (net in c(1,3)) {
        if (net == 1) { 
          name = 'PM'
        } else if (net == 3) { 
          name = 'AT' 
        }
        row = row+1
        meanConn$Subject[row] <- subjects[s]
        meanConn$Network[row] <- name
        meanConn$Region[row]  <- hipp
        # get mean of all seed and target connections between hippocampus and network (asymmetrical)
        meanConn$MeanConnectivity[row] <- mean(c(curMatrix[seeds == hipp,networks == net],
                                                 curMatrix[networks == net,seeds == hipp]))
      }
    }
  }
  meanConn$Subject <- as.factor(meanConn$Subject)
  meanConn$Network <- as.factor(meanConn$Network)
  meanConn$Region <- as.factor(meanConn$Region)
  
  return(meanConn)
}
# ------------------------------------------------------------------------------------------------------ #



## 3. Plot mean hipp to network connectivity change  ---------------------------------------------------- #
plot_hipp <- function(meanHipp){
  
  # plots retrieval vs encoding summary means:
  p <- ggplot(meanHipp, aes(x=Region, y=MeanConnectivity, fill=Network)) +
    stat_summary(fun.y = mean, geom="bar", alpha = 0.8, color = 'gray10', position=position_dodge(1)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6, alpha = 0.6, position=position_dodge(1)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96),
                 width = 0.45, color = "black", size = 0.65, position=position_dodge(1)) +
    scale_fill_manual(values = c('coral2','dodgerblue2')) + geom_hline(yintercept = 0) +
    ylab("Mean Beta") +
    ggtitle("Hipp Connectivity Strength") +
    theme(plot.title = element_text(hjust = 0.5, size=28), axis.line = element_line(colour = "black"), 
          axis.text = element_text(size = 22), axis.title = element_text(size = 24), panel.background = element_blank(),
          text = element_text(family="Helvetica")) 
  
  return(p)
}
# ------------------------------------------------------------------------------------------------ #



## 4. Calculate within/between network connectivity ---------------------------------------------- #
ppi_Network_connectivity <- function(subjects, ppiMatrix, networks){
  
  ### INPUT:
  # subjects  --> a list of subject IDs to include in analysis
  # ppiMatrix  --> a matrix of size seed x target x subject containing the psychologicalxphysiological beta estimates
  # networks   --> a vector of length seeds indicating network assignments

  ### OUTPUT:
  # network_ppi  --> the average connectivity estimate for network pairs per subject
  
  network_ppi = array(data = 0, c(length(networkNames),length(networkNames),length(subjects)))
  net <- unique(networks)
  
  for (sub in 1:length(subjects)) {
    curPPI = ppiMatrix[,,sub]
    diag(curPPI) <- NA
    
    for (seed in net) {
      for (targ in net) {
         netPPI = curPPI[networks == seed, networks == targ]
         network_ppi[seed,targ,sub]  <- mean(netPPI, na.rm = TRUE)
      } #end of loop through seed networks
    } #end of loop through target networks
  } #end of loop through subjects
  
  return(network_ppi)
}
# ------------------------------------------------------------------------------------------------ #



## 5. Get FDR corrected p-values per Network-Network Pair -------------------------------------------- #
get_network_pValues <- function(network_ppi, networkNames, nOrder, type){
  
  ### INPUT:
  # network_ppi  --> a matrix of size network x network x subject containing the mean psychologicalxphysiological beta estimates
  # networkNames --> a vector naming each of the unique networks (in index order)
  # nOrder    --> order to sort network levels
  # type      --> type of t-test - 'two-sided', 'greater', 'less'
  
  ### OUTPUT:
  # networkPvalues  --> a data.frame with the FDR corrected p value per network comparison and * to indicate < .05
  
  networkStats <- data.frame(array(NA, c(length(networkNames)*length(networkNames),9)))
  colnames(networkStats) <- c('seed','target','connectivity','SE','t','df','p','FDRp','sig')
  row = 0
  for (b in 1:nrow(network_ppi)) {
    for (a in 1:ncol(network_ppi)) {
      ppiVec <- network_ppi[a,b,]
      row = row + 1
      networkStats$seed[row] <- networkNames[a]
      networkStats$target[row] <- networkNames[b]
      networkStats$connectivity[row] <- mean(ppiVec)
      networkStats$SE[row] <- se(ppiVec)
      
      # one-sample t-test for this ROI-ROI connection across subjects
      stat <- t.test(ppiVec, mu=0, alternative = type)
      networkStats$t[row]  <- stat$statistic
      networkStats$df[row] <- stat$parameter
      networkStats$p[row]  <- stat$p.value
    } # end of loop through targets
  }  # end of loop through seeds

  networkStats$seed <- as.factor(networkStats$seed)
  networkStats$seed <- factor(networkStats$seed,levels(networkStats$seed)[nOrder])
  networkStats$target <- as.factor(networkStats$target)
  networkStats$target <- factor(networkStats$target,levels(networkStats$target)[nOrder])
  
  # adjust p values for matrix-level comparisons
  networkPAdjust <- p.adjust(networkStats$p, method = "fdr", n = length(networkStats$p))
  networkStats$FDRp <- networkPAdjust
  
  # add significance asterix for overlay on connection matrix (** FDR p < .05, * p < .05)
  for (r in 1:nrow(networkStats)) {
      if (as.numeric(networkStats$p[r]) < 0.05) {  #if significant, uncorrected, add single asterix
        networkStats$sig[r] <- '*'
      } 
      if (as.numeric(networkStats$FDRp[r]) < 0.05) {  #if FDR significant, add double asterix
        networkStats$sig[r] <- '**'
      }       
  }  # end of loop through seeds 
  
  return(networkStats)
}
# ------------------------------------------------------------------------------------------------ #



## 6. Plot mean network ppi matrix ----------------------------------------------------------------- #
plot_Network <- function(netStats, type){

  ### INPUT:
  # netStats --> means and p values for each seed-target comparison
  
  ### OUTPUT:
  # geom_tile ggplot
  
  myCol <- c("dodgerblue2","mediumorchid","firebrick2")
  
  p <- ggplot(data = netStats, aes(x=seed, y=target, fill=connectivity)) + 
    geom_tile(color = "white")+
    xlab('Seed') + ylab('Target') +
    geom_text(aes(label = sig, x=seed, y=target), size = 16, na.rm = TRUE, show.legend = FALSE) +
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=32, family="Helvetica", colour=myCol),
          axis.text.y = element_text(angle = 90, size=32, family="Helvetica", colour=myCol),
          axis.title = element_text(size=32,family="Helvetica"),
          panel.background = element_blank())
  
  if (type == 'greater') {
    limits = c(-0.4,0.8)
    p <- p + scale_fill_gradientn(colors = c("#6a51a3","#f7f7f7", "#ec7014"),
                                  limit = limits, space = "Lab", name="Change in\nConnectivity")
  } else {
    limits = c(-0.8,0.8)
    p <- p + scale_fill_gradient2(low = "#6a51a3", high = "#ec7014", mid = "#f7f7f7",
                                  midpoint = median(limits), limit = limits, space = "Lab", name="Change in\nConnectivity")
  }
  return(p)
}
# ------------------------------------------------------------------------------------------------ #



## 7. Get FDR corrected p-values per Seed-Target Pair -------------------------------------------- #
get_ppi_stats <- function(ppiMatrix, seeds, sOrder, type){

  ### INPUT:
  # ppiMatrix  --> a matrix of size seed x target x subject containing the psychologicalxphysiological beta estimates
  # seeds      --> a vector of seed names that were those used for the ppi analysis.
  # sOrder    --> order to sort ROI levels
  # type      --> type of t-test - 'two-tailed', 'greater', 'less'
  
  ### OUTPUT:
  # ppiStats  --> a data.frame with the FDR corrected p value per seed-target comparison and * to indicate < .05
  
  ppiPVals <- apply(ppiMatrix, c(1,2), function(x) t.test(x, mu = 0, alternative = type)$p.value)
  diag(ppiPVals) <- 1
  ppiTVals <- apply(ppiMatrix, c(1,2), function(x) t.test(x, mu = 0, alternative = type)$statistic)
  diag(ppiTVals) <- 0
  meanPPI <- apply(ppiMatrix, c(1,2), mean)
  
  ### melt p values
  ppiPVals <- data.frame(ppiPVals)
  colnames(ppiPVals) <- seeds
  ppiPVals$seeds <-seeds #seeds = rows
  ppiPVals$seeds <- as.factor(ppiPVals$seeds)
  ppiPVals$seeds = factor(ppiPVals$seeds,levels(ppiPVals$seeds)[sOrder])
  melted_ppiPVals <- melt(ppiPVals, id="seeds", na.rm = TRUE)
  
  ### melt t values
  ppiTVals <- data.frame(ppiTVals)
  colnames(ppiTVals) <- seeds
  ppiTVals$seeds <-seeds  #seeds = rows
  ppiTVals$seeds <- as.factor(ppiTVals$seeds)
  ppiTVals$seeds = factor(ppiTVals$seeds,levels(ppiTVals$seeds)[sOrder])
  melted_ppiTVals <- melt(ppiTVals, id="seeds", na.rm = TRUE)
  
  ## mean connectivity per connection:
  meanPPI <- data.frame(meanPPI)
  colnames(meanPPI) <- seeds
  meanPPI$seeds <-seeds  #seeds = rows
  meanPPI$seeds <- as.factor(meanPPI$seeds)
  meanPPI$seeds = factor(meanPPI$seeds,levels(meanPPI$seeds)[sOrder])
  melted_meanPPI <- melt(meanPPI, id="seeds", na.rm = TRUE)
  
  
  ### merge all connection stats ###
  melted_meanPPI <- merge(melted_meanPPI, melted_ppiTVals, by=c('seeds','variable'), sort=FALSE, suffixes = c('_conn','_t'))
  melted_meanPPI <- merge(melted_meanPPI, melted_ppiPVals, by=c('seeds','variable'), sort=FALSE)
  names <- colnames(melted_meanPPI)
  names[names == "value"] <- 'pvalue'
  names[names == 'variable'] <- 'targets'
  colnames(melted_meanPPI) <- names
  
  ## WHOLE MATRIX CORRECTION
  melted_meanPPI$pvalue[melted_meanPPI$seeds == melted_meanPPI$targets] <- NA #remove diagonal so not included in multiple comparison correction
  melted_ppiPAdjust <- p.adjust(melted_meanPPI$pvalue, method = "fdr", n = (sum(!is.na(melted_meanPPI$pvalue))))
  melted_meanPPI$pFDRAll <- melted_ppiPAdjust
  
  ### SEED LEVEL CORRECTION:
  for (s in 1:length(seeds)) {
    melted_ppiPAdjust <- p.adjust(melted_meanPPI$pvalue[melted_meanPPI$seeds == seeds[s]], method = "fdr", n = length(seeds)-1)
    melted_meanPPI$pFDRSeed[melted_meanPPI$seeds == seeds[s]] <- melted_ppiPAdjust #replace original p values with adjusted
  }
  
  # add significance asterix for overlay on connection matrix:
  melted_meanPPI$sig <- ''
  melted_meanPPI$sig[as.numeric(melted_meanPPI$pFDRAll) < 0.05] <- '*'
  
  return(melted_meanPPI)
}
# ------------------------------------------------------------------------------------------------ #



## 8. Plot mean ppi matrix -----------------------------------------------------------------------  #
plot_meanPPI <- function(ppiStats){
  
  ### INPUT:
  # ppiMatrix  --> a matrix of size seed x target x subject containing the psychologicalxphysiological beta estimates
  # melted_ppiPVals --> data frame of p values and * from get_FDR_pValues
  # seeds      --> a vector of seed names that were those used for the ppi analysis.
  # sOrder    --> order to sort ROI levels
  # title     --> title for plot

  ### OUTPUT:
  # geom_tile ggplot
  
  max <- round(max(abs(ppiStats$value_conn)),digits=3)
  min <- -(round(max(abs(ppiStats$value_conn)),digits=3))
  limits = c(min+(min*0.2),max+(max*0.2))
  
  myCol <- c("dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2",
             "mediumorchid","mediumorchid",
             "firebrick2","firebrick2","firebrick2","firebrick2","firebrick2")
  
  p <- ggplot(data = ppiStats, aes(x=seeds, y=targets, fill=value_conn)) + 
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("#6a51a3","#f7f7f7", "#ec7014"),
                         limit = limits, space = "Lab", name="Change in\nConnectivity") +
    xlab('Seed') + ylab('Target') +
    geom_text(aes(label = sig, x = seeds, y = targets), size = 12, na.rm = TRUE, show.legend = FALSE) +
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=20, colour=myCol), 
          axis.text.y = element_text(angle=45, vjust=0, hjust=1, size=20, colour=myCol), 
          axis.title = element_text(size=28,family="Helvetica"),
          panel.background = element_blank())
  
  return(p)
}
# ------------------------------------------------------------------------------------------------ #



## 9. plot invidual connection means -------------------------------------------- #
plot_seed_connections <- function(ppiMatrix, ppiStats, feature){
  
  # grabs subject data points for significant connections and plots:
  connections <- subset(ppiStats, pFDRSeed < .05 & (seeds == 'RSC' | seeds == 'PHC' | seeds == 'PRC' | seeds == 'AMYG'))
  nConn <- nrow(connections)  

  # get subject data per connection (only have mean value storedin ppiStats)
  nSub <- dim(ppiMatrix)[3]
  subData <- data.frame(array(NA, c(nSub*nConn,3)))
  colnames(subData) <- c('SubID','Connection','Beta')
  row = 0
  for (s in 1:nSub) {
    for (c in 1:nConn) {
      row = row+1
      seedIdx <- seeds == connections$seeds[c]
      targIdx <- seeds == connections$targets[c]
      subData$SubID[row] <- s
      subData$Connection[row] <- paste(connections$seeds[c],'-\n',connections$targets[c],sep = "")
      subData$Beta[row] <- ppiMatrix[seedIdx,targIdx,s]
    }
  }
  
  # plot:
  subData$Connection <- as.factor(subData$Connection)
  if (feature == 'Color') {
    subData$Connection <- factor(subData$Connection,levels(subData$Connection)[c(4,5,3,1,2)])
    myColor = 'coral2'
  } else if (feature == 'Scene'){
    subData$Connection <- factor(subData$Connection,levels(subData$Connection)[c(4,6,5,1,2,3)])
    myColor = 'dodgerblue2'
  }
  
  p <- ggplot(subData, aes(x=Connection, y=Beta)) +
    stat_summary(fun.y = mean, geom="bar", alpha = 0.5, fill=myColor, color=myColor, position=position_dodge(1)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, alpha = 0.8, fill=myColor, color=myColor, position=position_dodge(1)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), width = 0.45, color = "black", size = 0.65,
                 position=position_dodge(1)) +
    geom_hline(yintercept = 0) + ylab("Mean Beta") +
    ggtitle(paste('Connections Tracking\n',feature,'Precision',sep= " ")) +
    theme(plot.title = element_text(hjust = 0.5, size=28), axis.line = element_line(colour = "black"), 
          axis.text = element_text(size = 22), axis.title = element_text(size = 24), panel.background = element_blank(),
          legend.position = "none", text = element_text(family="Helvetica")) 
  
  return(p)
}

######################################################