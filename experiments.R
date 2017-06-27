library(wskm)
library(caret)
library(clue)
library(maxmatching)
library(SNFtool)
source('standard.R')
path <- "C:/Users/admin/Desktop/SZU/1_Dataset/Lymphoma.csv"

test_data <- read.csv(path, header = TRUE, sep = ',')



#######################################################################
#               Calculate groupinfo from data set                     #
#######################################################################
arr_to_group <- function(x){
  i <- 1
  feat_to_group <- c() 
  for (i in 1:length(x)){
    if (x[i] == 1){
      if (i %% feature_group_num == 0) 
        i <- feature_group_num
      else
        i <- i %% feature_group_num
      feat_to_group <- c(feat_to_group,i) 
      groupInfo <<- feat_to_group
      
    }
  }
  return (groupInfo)
}
#######################################################################
#                         match                                       #
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

#######################################################################
#                     Calculate benchmark                             #
#######################################################################
# parameter: alg <- 1 kmeans; 2 fgkmeans; 3  ewkmeans; 4 ## 
#           test_data
#           center: number of center
calculate_benchmark <- function(alg,test_data, center){
  
  rand_sum <- 0 
  accuracy_sum <- 0
  precision_sum <- 0
  recall_sum <- 0
  fmeasure_sum <- 0
  nmi_sum <- 0
  
  for(i in 1:5){
    cat('***************************************************************\n')
    cat('Calculating the', i,'time......\n')
    Sys.sleep(1)
    # seed=-1
    # if(seed<=0){
    #   seed <-runif(1,46,50)[1]
    # }
    # set.seed(42)
    if(alg == 1){
      cat('----------------------Using K-means algorithm--------------------------------\n')
      km <- kmeans(test_data, center)
    }
    else if(alg == 2){
      cat('----------------------Using FGK-means algorithm--------------------------------\n')
      groups <- arr_to_group()
      km <- fgkm(test_data,center)
    }
    else if(alg == 3){
      cat('----------------------Using EWK-means algorithm--------------------------------\n')
      km <- ewkm(test_data, center, maxrestart=-1)
    }
    else if(alg == 4){
      cat('----------------------Using EWK-means algorithm--------------------------------\n')
    }
    
    
    
    real.cluster <- c(rep(1,42), rep(2,9), rep (3,11))
    predict.cluster <- c(km$cluster)
    
    # matching clusters
    # minWeightBipartiteMatching(predict.cluster, real.cluster)
    #permuting predictive cluster
    
    #######################################################################
    #                matching and permuting cluster                       #
    #######################################################################
    matching <- minWeightBipartiteMatching(predict.cluster, real.cluster)
    clusterA <- predict.cluster  # map the labels from cluster A
    tmp <- sapply(1:length(matching), function(i) {
      clusterA[which(predict.cluster == i)] <<- matching[i]
    })
    clusterB <-  real.cluster
    
    
    cluster_table <- table(as.integer(clusterA), as.integer(clusterB))
  
    results <- confusionMatrix(cluster_table)
    
    #######################################################################
    #                            Rand Index                               #
    #######################################################################
    # std <- std.ext(as.integer(predict.cluster),as.integer(real.cluster))
    # rand_ind <- clv.Rand(std)
    # cat('---------------------------------------------------------------\n')
    # cat('Rand index is: ', rand_ind, '\n')
    # cat('---------------------------------------------------------------\n')
    rand_results <-com_accuracy(clusterA, clusterB,method = 4)
    rand_sum <- rand_sum + rand_results
    cat('The Rand Index between clusters ', ' is ',rand_results, '\n')
    cat('---------------------------------------------------------------\n')
    
    #######################################################################
    #                            accuracy                                 #
    #######################################################################
    overall.accuracy <- results$overall['Accuracy']
    accuracy_sum <- accuracy_sum + overall.accuracy
    cat('The overall accuracy is ', overall.accuracy, '\n')
    cat('---------------------------------------------------------------\n')
    #######################################################################
    #                  precision(AKA Pos Pred Value)                      #
    #######################################################################
    # precision_results <- c()
    # for(i in 1 : nrow(cluster_table)){
    #   precision_results[i] <- results$byClass[i, 5]
    #   cat('The precision of cluster ', i, ' is ',precision_results[i], '\n')
    # }
    precision_results <-com_accuracy(clusterA, clusterB,method = 1)
    precision_sum <- precision_sum + precision_results
    cat('The precision between clusters ', ' is ',precision_results, '\n')
    cat('---------------------------------------------------------------\n')
    
    #######################################################################
    #                   recall(AKA Sensitivity)                           #
    #######################################################################
    # recall_results <- c()
    # for(i in 1 : nrow(cluster_table)){
    #   recall_results[i] <- results$byClass[i, 1]
    #   cat('The recall of cluster ', i, ' is ', recall_results[i], '\n')
    # }
    recall_results <- com_accuracy(clusterA, clusterB,method = 2)
    recall_sum <- recall_sum + recall_results
    cat('The recall between two clusters', ' is ',recall_results, '\n')
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
    # f_measure <- c()
    # for(i in 1 : nrow(cluster_table)){
    #   f_measure[i] <- results$byClass[i, 7]
    #   cat('The F-measure of cluster ', i, ' is ', f_measure[i], '\n')
    # }
    f_measure <-  com_accuracy(clusterA, clusterB, mybeta = 1, method = 3)
    fmeasure_sum <- fmeasure_sum + f_measure
    cat('The f-measure between two clusters', ' is ',f_measure, '\n')
    cat('---------------------------------------------------------------\n')
    #######################################################################
    #              Normalized Mutual Information(NMI)                     #
    #######################################################################\
    # cat(predict.cluster,'\n', real.cluster, '\n', clusterA, '\n', clusterB)
    nmi <- calNMI(clusterA, clusterB)
    nmi_sum <- nmi_sum + nmi
    cat('The normalized mutual information is ',nmi,'\n')
    cat('---------------------------------------------------------------\n')
    Sys.sleep(1)
    cat('\n\n')
  }
  
  avg <- c(rand_sum/5, accuracy_sum/5, precision_sum/5, recall_sum/5, fmeasure_sum/5, nmi_sum/5)
  names(avg) <- c("Rand Index", "Accuracy", "Precision", "Recall", "F-measure", "NMI")
  print(avg)
  write.csv(avg, "C:/Users/admin/Desktop/Lymphoma_result.csv")
}

calculate_benchmark(1, test_data, 3)
