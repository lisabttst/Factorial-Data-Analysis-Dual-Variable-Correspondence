##### FONCTION QUI REALISE L'AFD #####

library(dplyr)


AFD <- function(data, Y_qual, Y_ind, diag_V=TRUE, norm=TRUE) {

  # n : nombre d'individus au total
  # p : nombre de variables quantitatives
  # q : nombre de variables qualitatives
  # nb_ind : nombre d'individus pour chaque variable qualitative
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  q <- length(Y_qual)
  nb_ind <- rep(0, q)
  
  
  # Calcule le nombre d'individus par catégorie
  for (i in 1:q) {
    nb_ind[i] <- nrow(data[data[,Y_ind] == Y_qual[i],])
  }
  
  
  # Centre les données
  for (col in 1:p) {
    data[,col] <- data[,col] - mean(data[,col])
  }
  
  
  # Réduit les données (oui par défaut)
  if(norm == TRUE) {
    for (i in 1:q) {
      for (col in 1:p) {
        data[data[,Y_ind] == Y_qual[i],][,col] <-  data[data[,Y_ind] == Y_qual[i],][,col] / sd(data[data[,Y_ind] == Y_qual[i],][,col])
      }
    }
  }
  
  
  # Moyennes sur chaque sous-population
  g_k <- matrix(0, q, p)
  for (row in 1:q) {
    qual_var <- Y_qual[row]
    for (col in 1:p) {
      g_k[row,][col] <- mean(data[data[,Y_ind] == qual_var,][,col]) 
    }
  }
  
  
  # Matrice B between, W within et V
  B <- matrix(0, p, p)
  for (i in 1:q) {
    B <- B + (nb_ind[i]/n) * g_k[i,] %*% t(g_k[i,])
  }
  V <- (1/n) * t(as.matrix(data[1:p])) %*% as.matrix(data[1:p])
  W <- V - B
  
  
  # Matrice C pour la diagonalisation
  C <- t(g_k)
  for (row in 1:nrow(C)) {
    for (col in 1:q) {
      C[row,col] <- sqrt(nb_ind[col]/n) * (C[row,col] - data[row,col]/n)
    }
  }
  print(C)

  
  
  # Diagonalise selon V ou W (V par défaut)
  if(diag_V == TRUE) {
    inv_V <- solve(V)
    l <- eigen(t(C) %*% inv_V %*% C)
    v <- inv_V %*% C %*% l$vectors
  } else {
    inv_V <- solve(W)
    l <- eigen(t(C) %*% inv_V %*% C)
    v <- inv_V %*% C %*% l$vectors
  }
  
  
  # Affiche les ratios de qualité de la projection
  print("== Qualité de la projection sur chaque axe ==")
  sum_vp <- sum(l$values)
  for (i in 1:q) {
    cat("Qualité sur l'axe", i, ":", l$values[i]/sum_vp, "\n")
  }
  
  print("== Contribution absolue des centres de gravités i à chacun des axes")
  for (i in 1:q) {
    for (axe in 1:(q-1)) {
      proportion <- nb_ind[i] / n
      cont <- proportion * (t(v[,axe]) %*% inv_V %*% g_k[i,])**2
      
      cat("Contribution absolue du groupe", i, "à l'axe", axe, ":",cont, "\n")
    }
  }
  
  print("== Contribution relative des centres de gravités i à chacun des axes")
  for (i in 1:q) {
    proportion <- nb_ind[i] / (n * l$values[i])
    cont <- proportion * (t(v[,i]) %*% inv_V %*% g_k[i,])**2
    
    cat("Contribution relative du groupe", i, "à l'axe", axe, ":",cont, "\n")
  }
  
  
  
  return(list("g_k"=g_k, "B"=B, "V"=V, "W"=W, "v"=v, "values"=l$values, "vectors"=l$vectors))
}
