################################################################################
# Homogeneity Tests for Related Samples ----                                                 
################################################################################
library(NADA2)
library(survival)

# Sign test ----                                                   
sign_related_samples <- function(df) {
  # Calculate differences
  differences <- df$x1 - df$x2
  non_zero_diffs <- differences[differences != 0]
  
  # Count positive and negative differences
  num_positive <- sum(non_zero_diffs > 0)
  num_negative <- sum(non_zero_diffs < 0)
  
  # Perform  test
  result <- binom.test(num_positive, num_positive + num_negative, p = 0.5)
  p_value <- result$p.value
  return(p_value)
}

## Paired Wilcoxon signed-rank test  ----                                                  
wilcoxon_signed_rank_related_samples <- function(df){
  # Adjust data
  df_paired <- data.frame(
    group1 = df$x1,
    group2 = df$x2
  )
  
  #Perform paired Wilcoxon signed-rank test
  result <- wilcox.test(df_paired$group1, df_paired$group2, paired = TRUE,
                        alternative = "two.sided")
  return(result$p.value)
}

## Akritas method ----                                                 
akritas_related_samples <- function(df) {
  # Fit survival curves for group 1 and group 2
  surv_km1 <- survfit(Surv(x1, delta1) ~ 1, data = df)
  surv_km2 <- survfit(Surv(x2, delta2) ~ 1, data = df)
  
  # Get survival probabilities and time
  surv1 <- surv_km1[["surv"]]
  surv2 <- surv_km2[["surv"]]
  time <- surv_km1[["time"]]
  
  # Store results in dataframe
  km_data <- data.frame(
    time = time,
    S1 = surv1,
    S2 = surv2,
    S_combined = (surv1 + surv2) / 2
  )
  
  # Calculate Akritas scores
  df <- df %>%
    rowwise() %>%
    mutate(
      S_group1 = approx(km_data$time, km_data$S1, x1, rule = 2)$y,
      S_group2 = approx(km_data$time, km_data$S2, x2, rule = 2)$y,
      S_combined = approx(km_data$time, km_data$S_combined, x1, rule = 2)$y,
      akritas_score1 = 1 - 0.5 * (S_group1 + delta1 * S_combined),
      akritas_score2 = 1 - 0.5 * (S_group2 + delta2 * S_combined)
    )
  
  # Calculate score diffrences
  df <- df %>%
    mutate(score_diff = akritas_score1 - akritas_score2)
  
  # Perform one-sample t-test for score differences
  result <- t.test(df$score_diff, mu = 0)
  return(result$p.value)
}


## PPW test ----                                              
ppw_related_samples <- function(df){
  # Perform data trasformations because test is designed for left-censored data
  df <- df %>% mutate(
    delta1 = ifelse(delta1 == 1, 0, 1),
    delta2 = ifelse(delta2 == 1, 0, 1)
  )
  
  max_value <- max(df$x1,df$x2)
  df_paired <- data.frame(
    group1 = max_value - df$x1,
    group2 = max_value - df$x2,
    delta1 = df$delta1,
    delta2 = df$delta2 
  )
  
  # Perform PPW test
  result <- ppw.test(df_paired$group1, df_paired$delta1,  df_paired$group2,
                     df_paired$delta2, alternative = "two.sided", printstat = FALSE)
  return(result$p.value)
}

## Weighted log-rank criterion (k=2) ----                                             
weighted_logrank_related_samples <- function(df){
  # Adjust data
  x <- as.matrix(df[, c("x1", "x2")])
  delta <- as.matrix(df[, c("delta1", "delta2")])
  # Sample size
  N <- nrow(x)
  # Weight function
  K <- 1/sqrt(N)
  
  ### Calculate V
  # initialize zero matrices y1 and y2
  y1 <- matrix(0, N, 2)
  y2 <- matrix(0, N, 2)
  
  x1_rep <- matrix(x[, 1], nrow = N, ncol = N, byrow = TRUE)
  x2_rep <- matrix(x[, 2], nrow = N, ncol = N, byrow = TRUE)
  
  y1[, 1] <- rowSums(x1_rep >= t(x1_rep))  
  y1[, 2] <- rowSums(x1_rep >= t(x2_rep))  
  y2[, 1] <- rowSums(x2_rep >= t(x1_rep))  
  y2[, 2] <- rowSums(x2_rep >= t(x2_rep))
  
  V1 <- K * y2[, 1] * delta[, 1] / (y1[, 1] + y2[, 1])
  V2 <- K * y1[, 2] * delta[, 2] / (y1[, 2] + y2[, 2])
  V <- (sum(V1) - sum(V2))
  
  ### Calculate C
  C1 <- (K^2) * y1[,1] * y2[,1] * delta[,1] / ((y1[,1] + y2[,1])^2)
  C2 <- (K^2) * y1[,2] * y2[,2] * delta[,2] / ((y1[,2] + y2[,2])^2)
  C <- (sum(C1) + sum(C2))
  
  ### Calculate R1 and R2
  R1 <- sqrt(N) * K * y2[,1] * delta[,1] / (y1[,1] + y2[,1])
  R2 <- sqrt(N) * K * y1[,2] * delta[,2] / (y1[,2] + y2[,2])
  
  ### Calculate R3 and R4
  y11 <- matrix(0, N, N)
  y12 <- matrix(0, N, N)
  y21 <- matrix(0, N, N)
  y22 <- matrix(0, N, N)
  
  for (l in 1:N) {
    y11[l, ] <- as.integer(x[, 1] >= x[l, 1]) 
    y12[l, ] <- as.integer(x[, 1] >= x[l, 2]) 
    y21[l, ] <- as.integer(x[, 2] >= x[l, 1]) 
    y22[l, ] <- as.integer(x[, 2] >= x[l, 2]) 
  }
  
  R3_1 <- sqrt(N) * K * delta[, 1] * y2[, 1] * y11 / ((y1[, 1] + y2[, 1])^2)
  R3_2 <- sqrt(N) * K * delta[, 2] * y2[, 2] * y12 / ((y1[, 2] + y2[, 2])^2)
  R3 <- colSums(R3_1) + colSums(R3_2)
  
  R4_1 <- sqrt(N) * K * delta[,1] * y1[,1] * y21 / ((y1[,1] + y2[,1])^2 )
  R4_2 <- sqrt(N) * K * delta[,2] * y1[,2] * y22 / ((y1[,2] + y2[,2])^2)
  R4 <- colSums(R4_1) + colSums(R4_2)
  
  ### Calculate C12
  C12 <- sum(R1 * R2 - R1 * R4 - R2 * R3 + R3 * R4) / N
  
  ### Calculate p-value
  sigma <- sqrt(C - 2 * C12) 
  Z_n <- V / sigma
  p_value <-  2 * (1 - pnorm(abs(Z_n)))
  return(p_value)
}

## Friedman test (k>2) ----                                              
friedman_related_samples <- function(df){
  # Adjust data
  df_paired <- data.frame(
    group1 = df$x1,
    group2 = df$x2,
    group3 = df$x3)
  
  # Perform Friedman test
  result <- friedman.test(as.matrix(df_paired))
  return(result$p.value)
}

## Weighted log-rank criterion  (k>2) ----                                            
weighted_logrank_related_samples_more_than_two <- function(df){
  # Adjust data
  x <- as.matrix(df[, c("x1", "x2", "x3")])  
  delta <- as.matrix(df[, c("delta1", "delta2", "delta3")])
  # Sample size
  N <- nrow(x)  
  # Number of samples
  k <- ncol(x)
  # Weight function
  K <- 1 / sqrt(N)  
  
  ### Calculate Y_i(X_lj) and Y(X_lj)
  # initialize zero matrices
  yi <- vector("list", length = k)
  for (g in 1:k) {
    yi[[g]] <- matrix(0, nrow = N, ncol = k)  #Y_i(X_lj)
  }
  y_total <- matrix(0, nrow = N, ncol = k) #Y(X_lj)
  
  for (j in 1:N) {  
    for (i in 1:N) {  
      for (g in 1:k) {  
        for (l in 1:k) {  
          if (x[i, g] >= x[j, l]) {
            yi[[g]][j, l] <- yi[[g]][j, l] + 1
          }
        }
      }
    }
    
    for (g in 1:k) {  
      for (f in 1:k) { 
        y_total[j, g] <- y_total[j, g] + yi[[f]][j, g]
      }
    }
    
  }
  
  ### Calculate V matrix
  # initialize zero matrix
  V <- matrix(0, nrow = 1, ncol = k - 1)
  
  for (g in 1:(k - 1)) {  
    for (j in 1:N) {  
      V_initial <- 0  
      for (l in 1:k) {  
        V_initial <- V_initial + K * delta[j, l] * yi[[g]][j, l] / y_total[j, l]
      }
      V[1, g] <- V[1, g] + delta[j, g] * K - V_initial
    }
  }
  
  #### Calculate C
  ## Calculate c_ij
  # initialize zero matrices
  sigma <- matrix(0, nrow = k - 1, ncol = k - 1)
  
  for (g in 1:(k - 1)) { 
    for (f in 1:(k - 1)) { 
      epsilon <- ifelse(g == f, 1, 0) # calculate epsilon_ij
      for (q in 1:N) {  
        for (s in 1:k) { 
          sigma[g, f] <- sigma[g, f] + K^2 * delta[q, s] * yi[[g]][q, s] / y_total[q, s] * (epsilon - yi[[f]][q, s] / y_total[q, s])
        }
      }
    }
  }
  
  ## Calculate c_ij(l,l')
  for (g in 1:(k - 1)) { 
    for (f in 1:(k - 1)) {  
      for (l in 1:k) {  
        epsilon_il <- ifelse(g == l, 1, 0)  # calculate epsilon_il
        
        for (l_llt in 1:k) { 
          epsilon_jl_llt <- ifelse(f == l_llt, 1, 0)  # calculate epsilon_jl_ll'
          
          if (l != l_llt) {  
            for (r in 1:N) {  
              
              r1_ilr <- delta[r, l] * sqrt(N) * K * (epsilon_il - yi[[g]][r, l] / y_total[r, l])
              r1_jl_llt_r <- delta[r, l_llt] * sqrt(N) * K * (epsilon_jl_llt - yi[[f]][r, l_llt] / y_total[r, l_llt])
              
              r2_ilr <- 0
              r2_jl_llt_r <- 0
              
              for (q in 1:N) {
                for (s in 1:k) {
                  if (x[r, l] >= x[q, s]) {
                    r2_ilr <- r2_ilr + delta[q, s] * sqrt(N) * K * (epsilon_il - yi[[g]][q, s] / y_total[q, s]) / y_total[q, s]
                  }
                  if (x[r, l_llt] >= x[q, s]) {
                    r2_jl_llt_r <- r2_jl_llt_r + delta[q, s] * sqrt(N) * K * (epsilon_jl_llt - yi[[f]][q, s] / y_total[q, s]) / y_total[q, s]
                  }
                }
              }
              sigma[g, f] <- sigma[g, f] + (1 / N) * (r1_ilr * r1_jl_llt_r - r1_ilr * r2_jl_llt_r - r1_jl_llt_r * r2_ilr + r2_ilr * r2_jl_llt_r)
            }
          }
        }
      }
    }
  }
  
  ### Calculate X2
  X2 <- V %*% solve(sigma) %*% t(V)
  
  ### Calculate p-value (the upper tail probability)
  p_value <- pchisq(X2, df = k - 1, lower.tail = FALSE) 
  return(p_value)
}
