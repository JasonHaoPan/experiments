com_accuracy <- function(vec1, vec2, method = 0, mybeta = 0) {
  
  ## vec1: the cluster result by experiment
  ## vec2: the real cluster, stand cluster
  ## method=0 purity
  ## method=1 precision
  ## method=2 recall
  ## method=3 F-score with para mybeta
  ## method=4 RI
  
  ## modify-time 2016-1-3 fix the bug when there is only one data point in a class in vec1 or vec2
  
  if (class(vec1) == "communities") {
    vec1 = membership(vec1)
    
  }
  
  if (class(vec2) == "communities") {
    vec2 = membership(vec2)
    
  }
  
  n = length(vec1)
  
  nodes = 1:n
  
  K1 = length(unique(vec1))
  K2 = length(unique(vec2))
  com1 = c()
  com2 = c()
  for (i in 1:K1) {
    com1[i] = list(nodes[vec1 == i])
  }
  for (i in 1:K2) {
    com2[i] = list(nodes[vec2 == i])
  }
  if (method == 0) {
    ## purity
    cor_num = 0
    may_id = 1:K2
    for (i in 1:K1) {
      
      temp = sapply(may_id, function(x) {
        return(length(intersect(com1[[i]], com2[[x]])))
        
      })
      cor_num = cor_num + max(temp)
      
    }
    
    accuracy = cor_num/n
    return(accuracy)
  }
  if (method != 0) 
  {
    ## precision
    expe_cluster_mat = matrix(0, n, n)
    real_cluster_mat = matrix(0, n, n)
    for (i in 1:K1) {
      # if there is one node in the cluster, next
      if (length(com1[[i]]) == 1) {
        next
      }
      nodepair_id = combn(com1[[i]], 2)
      node_id = apply(nodepair_id, 2, function(x) {
        x = sort(x)
        x[1] + (x[2] - 1) * n
      })
      expe_cluster_mat[node_id] = 1
    }
    expe_cluster_mat[lower.tri(expe_cluster_mat)] = 2
    for (i in 1:K2) {
      # if there is one node in the cluster, next for stand cluster, there should not happen
      if (length(com2[[i]]) == 1) {
        next
      }
      nodepair_id = combn(com2[[i]], 2)
      node_id = apply(nodepair_id, 2, function(x) {
        x = sort(x)
        x[1] + (x[2] - 1) * n
      })
      real_cluster_mat[node_id] = 1
    }
    real_cluster_mat[lower.tri(real_cluster_mat)] = 2
    diag(expe_cluster_mat) = 2
    diag(real_cluster_mat) = 2
    exp_p = which(expe_cluster_mat == 1)
    rel_p = which(real_cluster_mat == 1)
    exp_n = which(expe_cluster_mat == 0)
    rel_n = which(real_cluster_mat == 0)
    TP = length(intersect(rel_p, exp_p))
    FP = length(intersect(rel_n, exp_p))
    
    ## -n, delete the effect by the diag element, which both be 0
    TN = length(intersect(rel_n, exp_n))
    FN = length(intersect(rel_p, exp_n))
    P = TP/(TP + FP)
    R = TP/(TP + FN)
  }  ## end if method !=0
  
  
  # cat(TP, TN, FP, FN, "\n")
  
  if (method == 1) {
    accuracy = P
  } else if (method == 2) {
    accuracy = R
  } else if (method == 3) {
    accuracy = (mybeta + 1) * P * R/(mybeta * P + R)
  } else if (method == 4) {
    accuracy = (TP + TN)/(TP + TN + FP + FN)
  }
  return(accuracy)
}