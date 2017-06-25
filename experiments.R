library(wskm)
library(caret)
library(clue)
library(maxmatching)
path <- "C:/Users/admin/Desktop/SZU/1_Dataset/Lymphoma.csv"

test_data <- read.csv(path, header = TRUE, sep = ',')

#######################################################################
#                         matching                                    #
#######################################################################
# labels from cluster A will be matched on the labels from cluster B
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
  require(clue)
  idsA <- unique(clusteringA)  # distinct cluster ids in a
  idsB <- unique(clusteringB)  # distinct cluster ids in b
  nA <- length(clusteringA)  # number of instances in a
  nB <- length(clusteringB)  # number of instances in b
  if (length(idsA) != length(idsB) || nA != nB) {
    stop("number of cluster or number of instances do not match")
  }
  
  nC <- length(idsA)
  tupel <- c(1:nA)
  
  # computing the distance matrix
  assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
  for (i in 1:nC) {
    tupelClusterI <- tupel[clusteringA == i]
    solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
      nA_I <- length(tupelA_I)  # number of elements in cluster I
      tupelB_I <- tupel[clusterIDsB == i]
      nB_I <- length(tupelB_I)
      nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
      return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
    }, clusteringB, tupelClusterI)
    assignmentMatrix[i, ] <- solRowI
  }
  
  # optimization
  result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
  attr(result, "assignmentMatrix") <- assignmentMatrix
  return(result)
}
sum <- 0;
for(i in 1:5){
  cat('***************************************************************\n')
  cat('Calculating the', i,'time......\n')
  Sys.sleep(2)
  seed=-1
  if(seed<=0){
    seed <-runif(1,0,10000000)[1]
  }
  set.seed(seed)
  km <- ewkm(test_data, 3, maxrestart=-1)
  real.cluster <- c(rep(1,42), rep(2,9), rep (3,11))
  #real.cluster <- test_data
  predict.cluster <- c(km$cluster)
  # calculate rand index
  std <- std.ext(as.integer(predict.cluster),as.integer(real.cluster))
  rand_ind <- clv.Rand(std)
  cat('---------------------------------------------------------------\n')
  cat('Rand index is: ', rand_ind, '\n')
  cat('---------------------------------------------------------------\n')
  sum <- sum+rand_ind
  
  # matching clusters
  # minWeightBipartiteMatching(predict.cluster, real.cluster)
  #permuting predictive cluster
  
  #######################################################################
  #                matching and permuting cluster                       #
  #######################################################################
  matching <- minWeightBipartiteMatching(predict.cluster, real.cluster)
  matching
  clusterA <- predict.cluster  # map the labels from cluster A
  tmp <- sapply(1:length(matching), function(i) {
    clusterA[which(predict.cluster == i)] <<- matching[i]
  })
  clusterB <-  real.cluster
  
  # clusterA
  # clusterB
  
  cluster_table <- table(as.integer(clusterA), as.integer(clusterB))
  cluster_table
  results <- confusionMatrix(cluster_table)
  #######################################################################
  #                            accuracy                                 #
  #######################################################################
  overall.accuracy <- results$overall['Accuracy']
  cat('The overall accuracy is ', overall.accuracy, '\n')
  cat('---------------------------------------------------------------\n')
  #######################################################################
  #                  precision(AKA Pos Pred Value)                      #
  #######################################################################
  for(i in 1 : nrow(cluster_table)){
    precision[i] <- results$byClass[i, 5]
    cat('The precision of cluster ', i, ' is ', precision[i], '\n')
  }
  cat('---------------------------------------------------------------\n')
  #######################################################################
  #                   recall(AKA Sensitivity)                           #
  #######################################################################
  for(i in 1 : nrow(cluster_table)){
    recall[i] <- results$byClass[i, 1]
    cat('The recall of cluster ', i, ' is ', recall[i], '\n')
  }
  cat('---------------------------------------------------------------\n')
  #######################################################################
  #                     F-measure(AKA  F1 Score)                        #
  #######################################################################
  
  #use my own function
  # for(i in 1 : nrow(cluster_table))
  # {
  #   f_measure[i] <- 0;
  #   f_measure[i] <- f_measure[i] + 2*precision[i]*recall[i]/(precision[i]+recall[i])
  #   cat('The F-measure of cluster ', i, ' is ', f_measure[i], '\n')
  # }
  
  for(i in 1 : nrow(cluster_table)){
    f_measure[i] <- results$byClass[i, 7]
    cat('The F-measure of cluster ', i, ' is ', f_measure[i], '\n')
  }
  cat('---------------------------------------------------------------\n')
  Sys.sleep(2)
  cat('\n\n')
}
avg <- sum/5
print(avg)


