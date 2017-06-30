library(wskm)
library(caret)
library(clue)
library(maxmatching)
library(SNFtool)
library(kernlab)
library(speccalt)
source('standard.R')
source('matching.R')
# path <- "C:/Users/admin/Desktop/SZU/1_Dataset/Lymphoma.csv"



#######################################################################
#                      Initialize group                              #
#######################################################################
initialize_group <- function(feature_num,group_num){
  my_group <- sample(1:group_num,feature_num,replace = TRUE)
}


#######################################################################
#                     Calculate benchmark                             #
#######################################################################
# parameter: alg <- 1 kmeans; 2 fgkmeans; 3  ewkmeans; 4 ## 
#           test_data
#           center: number of center
calculate_benchmark <- function(alg,test_data, center, real_cluster){
  
  rand_sum <- 0 
  accuracy_sum <- 0
  precision_sum <- 0
  recall_sum <- 0
  fmeasure_sum <- 0
  nmi_sum <- 0
  
  for(i in 1:5){
    # set.seed(42)
    if(alg == 1){
      # cat('----------------------Using K-means algorithm--------------------------------\n')
      km <- kmeans(test_data, center)
      predict.cluster <- c(km$cluster)
    }else if(alg == 2){
      # cat('----------------------Using FGK-means algorithm--------------------------------\n')
      grouping <- initialize_group(ncol(test_data),3)
      km <- fgkm(test_data, center, grouping, 1, 1)
      predict.cluster <- c(km$cluster)
    }else if(alg == 3){
      # cat('----------------------Using EWK-means algorithm--------------------------------\n')
      km <- ewkm(test_data, center, maxrestart=-1)
      predict.cluster <- c(km$cluster)
    }else if(alg == 4){
      # cat('----------------------Using TWK-means algorithm--------------------------------\n')
      grouping <- initialize_group(ncol(test_data), 3)
      km <- twkm(test_data, center,grouping, 1, 1)
      predict.cluster <- c(km$cluster)
    }else if(alg == 5){
      # cat('----------------------Using Spectral Clustering algorithm--------------------------------\n')
      km <- speccl(test_data, center)
      predict.cluster <- c(km$clusters)
    }
    cat('***************************************************************\n')
    cat('Calculating the', i,'time......\n')
    Sys.sleep(1)
    real.cluster <- real_cluster
   
    # matching clusters
    minWeightBipartiteMatching(predict.cluster, real.cluster)
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
    accuracy_results <-com_accuracy(clusterA, clusterB,method = 0)
    accuracy_sum <- accuracy_sum + accuracy_results
    cat('The overall accuracy is ', accuracy_results, '\n')
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
    f_measure <-  com_accuracy(clusterA, clusterB, mybeta = 1, method = 3)
    fmeasure_sum <- fmeasure_sum + f_measure
    cat('The f-measure between two clusters', ' is ',f_measure, '\n')
    cat('---------------------------------------------------------------\n')
    
    #######################################################################
    #              Normalized Mutual Information(NMI)                     #
    #######################################################################\
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
  Sys.sleep(3)
  
  return(avg)
}

filenames <- list.files("C:/Users/admin/Desktop/SZU/1_Dataset/", pattern = ".csv")
for(i in 1:length(filenames))
{
  path <- paste("C:/Users/admin/Desktop/SZU/1_Dataset/",filenames[i], sep = "")
  test_data <- read.csv(as.character(path), header = TRUE, sep = ',')
  cat(as.character(path),"\n")
  Sys.sleep(3)
  if(filenames[i]=="CNS.csv"||
      filenames[i]=="Colon.csv"||
        filenames[i]=="Leukemia_3c.csv"||
          filenames[i]=="Leukemia_4c.csv"){
    
    col <- ncol(test_data)
    if(filenames[i]=="CNS.csv"){
      real_cluster <- as.numeric(test_data[,col]+1)
    }else{
      real_cluster <- as.numeric(test_data[,col])
    }
    
    test_data <- test_data[,-col]
  }else{
    if(filenames[i]=="LeukemiaTXT.csv"){
      real_cluster <- as.numeric(test_data[,1])
    }else{
      real_cluster <- (test_data[,1]+1)
    }
    test_data <- test_data[,-1]
  }
  table(real_cluster)
  center <- length(table(real_cluster))
  bench_result <- c()
  for(j in 1:5)
  {
    if(j == 1)
    {
      cat("******************Using K-means Algorithm*********************\n")
      kmeans_alg <- calculate_benchmark(j, test_data, center, real_cluster)
      bench_result <- cbind(bench_result,kmeans_alg)
    }else if(j == 2)
    {
      cat("*****************Using FGK-means Algorithm*********************\n")
      # fgkmeans_alg <- calculate_benchmark(j, test_data, center, real_cluster)
      # bench_result <- cbind(bench_result,fgkmeans_alg)
    }else if(j == 3)
    {
      cat("******************Using EWK-means Algorithm********************\n")
      ewkmeans_alg <- calculate_benchmark(j, test_data, center, real_cluster)
      bench_result <- cbind(bench_result,ewkmeans_alg)
    }else if(j == 4)
    {
      cat("******************Using TWK-means Algorithm*********************\n")
      twkmeans_alg <- calculate_benchmark(j, test_data, center, real_cluster)
      bench_result <- cbind(bench_result,twkmeans_alg)
    }
    else if(j == 5)
    {
      cat("***************Using Spectral Clustering Algorithm***************\n")
      spectral_alg <- calculate_benchmark(j, test_data, center, real_cluster)
      bench_result <- cbind(bench_result,spectral_alg)
    }
  }
  outputpath <- paste("C:/Users/admin/Desktop/result/",filenames[i],sep = "")
  cat("Writing to file...\n")
  Sys.sleep(3)
  write.csv(bench_result, outputpath)
  cat("Done\n")
}



# library(clusterSim)
# 
# test <- read.csv("C:/Users/admin/Desktop/SZU/1_Dataset/adenocarcinoma.csv", header = TRUE, sep = ",")
# test <- test[, -1]
# 
# res <- speccl(test, nc = 2)
# res$clusters

