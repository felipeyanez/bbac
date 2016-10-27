###############################################################################
### BREGMAN BLOCK AVERAGE CO-CLUSTERING ALGORITHM (BANERJEE ET AL., 2007)
###############################################################################

bbac <- function(Z, k, l, W = NULL, distance = "euclidean", scheme = 6, 
                 errobj = 1e-6, niters = 100, nruns = 5, epsilon = 1e-8) {
  
  error <- Inf
  error_now <- Inf
  dist <- pmatch(tolower(distance), c("euclidean","divergence")) - 1
  
  for (r in 1:nruns) {
    
    # Initialization of R and C
    R <- diag(k)[sample(k, dim(Z)[1], replace = TRUE),]
    C <- diag(l)[sample(l, dim(Z)[2], replace = TRUE),]
    
    for (s in 1:niters) {
      
      # Row estimation
      rs <- get_scheme(dist, scheme, "row", R, Z, C, W, epsilon)
      ra <- assign_cluster(dist,    Z,    rs$Zrowc,    rs$Zrowv, W, epsilon)
      R  <- ra$Cluster
      
      # Column estimation
      cs <- get_scheme(dist, scheme, "col", R, Z, C, W, epsilon)
      ca <- assign_cluster(dist, t(Z), t(cs$Zcolc), t(cs$Zcolv), W, epsilon)
      C  <- ca$Cluster
      
      # 
      if (abs(ca$Error - error_now) < errobj) {
        status <- paste("converged in",s,"iterations")
        return(list(R = R, C = C, status = status))
      }
      
      # Update objective value
      error_now <- ca$Error
      
    }
    
    # Keep pair with min error
    if (error_now < error) {
      R_star <- R
      C_star <- C
      error <- error_now
    }
    
  }
  
  status <- paste("reached maximum of",niters,"iterations")
  return(list(R = R_star, C = C_star, status = status))
  
}

###############################################################################
### HELPER FUNCTIONS 
###############################################################################

get_scheme <- function(dist, scheme, row_col, R, Z, C, W, epsilon) {
  if (scheme==1) {
    return(scheme1(dist, row_col, R, Z, C, W, epsilon))
  } else if (scheme==2) {
    return(scheme2(dist, row_col, R, Z, C, W, epsilon))
  } else if (scheme==3) {
    return(scheme3(dist, row_col, R, Z, C, W, epsilon))
  } else if (scheme==4) {
    return(scheme4(dist, row_col, R, Z, C, W, epsilon))
  } else if (scheme==5) {
    return(scheme5(dist, row_col, R, Z, C, W, epsilon))
  } else if (scheme==6) {
    return(scheme6(dist, row_col, R, Z, C, W, epsilon))
  }
}

scheme1 <- function(dist, row_col, R, Z, C, W, epsilon) {
  e1 <- matrix(1, dim(Z)[1])
  e2 <- matrix(1, dim(Z)[2])
  RCavg <- calculate_average( R, Z, e2, W, epsilon)
  CCavg <- calculate_average(e1, Z,  C, W, epsilon)
  if (row_col=="row") {
    A <- c("(e1 %*% CCavg %*% t(C))", "(array(mean(Z), dim(Z)))")
    return(list(Zrowc = op1(dist, A), Zrowv = RCavg %*% t(e2)))
  } else if (row_col=="col") {
    A <- c("(R %*% RCavg %*% t(e2))", "(array(mean(Z), dim(Z)))")
    return(list(Zcolc = op1(dist, A), Zcolv = e1 %*% CCavg))
  }
}

scheme2 <- function(dist, row_col, R, Z, C, W, epsilon) {
  CoCavg <- calculate_average(R, Z, C, W, epsilon)
  if (row_col=="row") {
    return(list(Zrowc = array(dist, dim(Z)), Zrowv = CoCavg %*% t(C)))
  } else if (row_col=="col") {
    return(list(Zcolc = array(dist, dim(Z)), Zcolv = R %*% CoCavg))
  }
}

scheme3 <- function(dist, row_col, R, Z, C, W, epsilon) {
  I1 <- diag(dim(Z)[1])
  e2 <- matrix(1, dim(Z)[2])
  Ravg   <- calculate_average(I1, Z, e2, W, epsilon)
  RCavg  <- calculate_average( R, Z, e2, W, epsilon)
  CoCavg <- calculate_average( R, Z,  C, W, epsilon)
  if (row_col=="row") {
    A <- c("(CoCavg %*% t(C))", "(RCavg %*% t(e2))")
    return(list(Zrowc = Ravg %*% t(e2), Zrowv = op1(dist, A)))
  } else if (row_col=="col") {
    A <- c("(Ravg %*% t(e2))", "(R %*% RCavg %*% t(e2))")
    return(list(Zcolc = op1(dist, A), Zcolv = R %*% CoCavg))
  }
}

scheme4 <- function(dist, row_col, R, Z, C, W, epsilon) {
  I2 <- diag(dim(Z)[2])
  e1 <- matrix(1, dim(Z)[1])
  Cavg   <- calculate_average(e1, Z, I2, W, epsilon)
  CCavg  <- calculate_average(e1, Z,  C, W, epsilon)
  CoCavg <- calculate_average( R, Z,  C, W, epsilon)
  if (row_col=="row") {
    A <- c("(e1 %*% Cavg)", "(e1 %*% CCavg %*% t(C))")
    return(list(Zrowc = op1(dist, A), Zrowv = CoCavg %*% t(C)))
  } else if (row_col=="col") {
    A <- c("(R %*% CoCavg)", "(e1 %*% CCavg)")
    return(list(Zcolc = e1 %*% Cavg, Zcolv = op1(dist, A)))
  }
}

scheme5 <- function(dist, row_col, R, Z, C, W, epsilon) {
  I1 <- diag(dim(Z)[1])
  I2 <- diag(dim(Z)[2])
  e1 <- matrix(1, dim(Z)[1])
  e2 <- matrix(1, dim(Z)[2])
  Ravg   <- calculate_average(I1, Z, e2, W, epsilon)
  Cavg   <- calculate_average(e1, Z, I2, W, epsilon)
  RCavg  <- calculate_average( R, Z, e2, W, epsilon)
  CCavg  <- calculate_average(e1, Z,  C, W, epsilon)
  CoCavg <- calculate_average( R, Z,  C, W, epsilon)
  A <- "op2(dist, c((Ravg %*% t(e2)), (e1 %*% Cavg)))"
  if (row_col=="row") {
    A <- c(A, "(e1 %*% CCavg %*% t(C))")
    B <- c("(CoCavg %*% t(C))", "(RCavg %*% t(e2))")
    return(list(Zrowc = op1(dist, A), Zrowv = op1(dist, B)))
  } else if (row_col=="col") {
    A <- c(A, "(R %*% RCavg %*% t(e2))")
    B <- c("(R %*% CoCavg)", "(e1 %*% CCavg)")
    return(list(Zcolc = op1(dist, A), Zcolv = op1(dist, B)))
  }
}

scheme6 <- function(dist, row_col, R, Z, C, W, epsilon) {
  I1 <- diag(dim(Z)[1])
  I2 <- diag(dim(Z)[2])
  CoCavg  <- calculate_average( R, Z,  C, W, epsilon)
  RC_Cavg <- calculate_average( R, Z, I2, W, epsilon)
  CC_Ravg <- calculate_average(I1, Z,  C, W, epsilon)
  if (row_col=="row") {
    A <- c("(RC_Cavg)", "(CoCavg %*% t(C))")
    return(list(Zrowc = CC_Ravg %*% t(C), Zrowv = op1(dist, A)))
  } else if (row_col=="col") {
    A <- c("(CC_Ravg)", "(R %*% CoCavg)")
    return(list(Zcolc = R %*% RC_Cavg, Zcolv = op1(dist, A)))
  }
}

calculate_average <- function(Left, Z, Right, W, epsilon) {
  if (is.null(W)) {W <- array(1, dim(Z))} else {Z <- W * Z}
  numerator <- t(Left) %*% Z %*% Right + mean(Z) * epsilon
  denominator <- t(Left) %*% W %*% Right + epsilon
  return(numerator/denominator)
}

op1 <- function(dist, A) {
  B <- paste(A[1], ifelse(dist==0, "-", "/"), A[2])
  return(eval.parent(parse(text = B), n = 1))
}

op2 <- function(dist, A) {
  B <- paste(A[1], ifelse(dist==0, "+", "*"), A[2])
  return(eval.parent(parse(text = B), n = 1))
}

assign_cluster <- function(dist, Z, X, Y, W, epsilon) {
  D <- similarity_measure(dist, Z, X, Y, W, epsilon)
  id <- sapply(1:dim(D)[1], function(i) sort(D[i,], index.return = TRUE)$ix[1])
  res <- sapply(1:dim(D)[1], function(i) sort(D[i,])[1]^(2-dist))
  return(list(Cluster = diag(dim(Y)[1])[id,], Error = sum(res)))
}

similarity_measure <- function(dist, Z, X, Y, W, epsilon) {
  if (is.null(W)) W <- array(1, dim(Z))
  if (dist==0) {
    euc <- function(i) rowSums(W * (Z - X - rep(Y[i,], each = dim(Z)[1]))^2)
    return(sapply(1:dim(Y)[1], euc))
  } else if (dist==1) {
    return((W * X) %*% t(Y + epsilon) - (W * Z) %*% log(t(Y + epsilon)))
  }
}

plot_coclusters <- function(Z, R, C) {
  # Sort matrix
  Y <- t(Z[(R * (1:nrow(R)))[R != 0], (C * (1:nrow(C)))[C != 0]])
  # Plot sorted matrix
  image(seq(0, 1, length.out = dim(Y)[1]), seq(0, 1, length.out = dim(Y)[2]),
        Y, col = grey((0:12)/12), axes = FALSE, xlab = "", ylab = "")
  # Print row clusters
  row_clust <- (head(cumsum(colSums(R)), -1) - 0.5)/(ncol(Y) - 1)
  invisible(sapply(1:length(row_clust), function(i) 
    segments(-0.5, row_clust[i], 1.5, row_clust[i], col = 2, lwd = 2)))
  # Print column clusters
  col_clust <- (head(cumsum(colSums(C)), -1) - 0.5)/(nrow(Y) - 1)
  invisible(sapply(1:length(col_clust), function(i) 
    segments(col_clust[i], -0.5, col_clust[i], 1.5, col = 2, lwd = 2)))
}