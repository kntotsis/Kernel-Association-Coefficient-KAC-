# Kernel Association Coefficient (KAC) function
# x = a numeric vector, matrix or dataframe;
## if ncol(x) > 1, then x must be a dataframe.

# y = a numeric vector, matrix or dataframe (optional);
## if y is missing, then x must be a dataframe and ncol(x) >= 2.
## if y is not missing then ncol(y) must be 1.

# const = a numeric single constant value for data readjustment (optional);
## if x (or y or both) contain zero values then:
### if const is missing, then x = x + 10^(-3)
### if const is not missing, then x = x +10^(-const)


KAC <- function(x, y, const) {
  options(warn = -1) #suppress the warnings
  
  if (missing(y)) {
    # if y is missing
    if (is.data.frame(x) == FALSE) {
      # and if x is not a dataframe, to return:
      stop(paste("input data must be a numeric dataframe"))
    } else if (ncol(x) < 2) {
      # if x is indeed a dataframe but without at least two features, to return:
      stop(paste("x must contain at least two columns when y is missing"))
    } else {
      # if x contains at least two features, rename is as `data`
      data <- x
    }
  } else {
    # if y is not missing
    if (length(y) == 0) {
      # it can't be an empty vector
      stop("y must not be empty")
    } else if (ncol(y) != 1) {
      # can't have more than one column
      stop("y can not have more than one columns")
    } else{
      # combine it with x in the dataframe `data`
      data <- as.data.frame(cbind(x, y))
    }
  }
  
  if (any(data == 0, na.rm = TRUE)) {
    # if `data` contains zeros then:
    if (missing(const)) {
      # if the optional argument `const` is omitted, then add 10^(-3) in all elements of `data`, otherwise
      data <- data + 10 ^ (-3)
    } else {
      if (!is.numeric(const) || length(const) > 1) {
        stop(paste("const must be a single numeric value"))
      }
      # if the optional argument `const` is given, then add 10^(-const) in all elements of `data`.
      data <- data + 10 ^ (-const)
    }
  }
  
  
  # Step 1: calculation of the ceontering matrix C_m (Eq. 2)
  m <- nrow(data) # number of observatioåns
  I <- diag(m) # m-size identity matrix
  one_m <-
    as.matrix(rep(c(1), times = m)) # 1m = (m x 1): vector of ones
  
  C_m <- I - (1 / m) * one_m %*% t(one_m)
  
  # Step 2: Compute the uncentered Grammian matrices (G^Z.) (Eq. 4)
  G <-
    lapply(1:ncol(data), function(i)
      ((1 / var(data[, i])) * abs(outer(data[, i], data[, i]))) ^ 2)
  
  # Step 3: Compute the centered Grammian matrices (G_c^Z.) (Eq. 5)
  G_c <- lapply(1:ncol(data), function(i)
    C_m * G[[i]] * C_m)
  
  # Step 4a: Compute the w. quantity (Eq. 7)
  w. <- t(matrix(sapply(G_c, \(x) min(colSums(
    x
  )))))
  
  # Step 4b: Compute the pivot matrices (P..) (Eq. 6)
  P_i.j <- # P_i.j is the pivot matrixes when i≠j
    combn(
      seq_along(G_c),
      2,
      FUN = function(i)
        
        (-0.5 * log(abs(
          I - (min(w.[, i]) ^ (-2)) *
            G_c[[i[1]]] * G_c[[i[2]]]
        ))),
      simplify = FALSE
    )
  names(P_i.j) <-
    paste0("P(", combn(seq_along(G_c), 2, FUN = paste, collapse = ","), ")")
  
  # Step 5: Compute the Kernel Association Matrix (Eq. 11)
  P_i.i <- # P_i.i is the pivot matrixes when i=j
    lapply(seq_along(G_c), \(i)  (-0.5 * log(abs(
      I - (w.[, i] ^ (-2)) * (G_c[[i]] * G_c[[i]])
    ))))
  names(P_i.i) <-
    paste0("P(", seq(1, ncol(data)), ",", seq(1, ncol(data)), ")")
  
  P.. <- list() # P.. is the complete list of pivot matrices
  for (i in 1:ncol(data)) {
    P..[[paste0("P(", i, ",", i, ")")]] <-
      P_i.i[[paste0("P(", i, ",", i, ")")]]
    for (j in 1:ncol(data)) {
      if (i != j) {
        P..[[paste0("P(", i, ",", j, ")")]] <-
          P_i.j[[paste0("P(", i, ",", j, ")")]]
      }
    }
  }
  
  P_lambdas <- do.call(cbind, lapply(P.., function(x)
    round(eigen(x)$values, 2))) # the λ values
  
  P_lambdas <-
    matrix(P_lambdas[1,],
           ncol = ncol(P_lambdas),
           dimnames = list(NULL, colnames(P_lambdas))) # the highest (first) λ value
  
  KAC <- matrix(rep(0, ncol(data) ^ 2), nrow = ncol(data))
  KAC[lower.tri(KAC, diag = TRUE)] <- P_lambdas
  KAC[upper.tri(KAC)] <- t(KAC)[upper.tri(KAC)] # the unprojected KAC matrix (Eq. 11) 
  
  # Step 6: Compute the projected Kernel Association Matrix (Eq. 13)
  KAC <- KAC /
    (sqrt(diag(KAC)) %*% t(sqrt(diag(KAC)))) #Eq. 12
  dimnames(KAC) <- list(colnames(data), colnames(data))
  
  # Return the KAC matrix
  if (ncol(data) == 2) {
    KAC <- KAC[1, 2] # the association between thw 2 features
  }
  if (any(KAC < 0)) {
    KAC[KAC < 0] <- 0
    message(
      "Values close to zero but slightly negative have been detected in the KAC matrix and have been set to 0. This occurrence is likely a result of inadequate data preparation, leading to computational errors. It is recommended to thoroughly review and preprocess the input data for more efficient results."
    )
    return(KAC)
  } else {
    return(KAC)
  }
}
