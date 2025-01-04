################################################################################
# Simulated examples for testing homogeneity of related samples ----
################################################################################
library(dplyr)
library(FAdist)
library(copula)

# Complete data, two-sample case ----
sample_size <- 500 # sample size
tau_target <- 2/3 # Kendall tau
lambda1 <- 1 # scale parameter for t1
lambda2 <-1 # scale parameter for t2
nu1 <- 1 # shape parameter for t1
nu2 <- 1 # shape parameter for t2
theta1 <- 0 # location parameter for t1
theta2 <- 0 # location parameter for t2 

# Generate data from the copula (for other copulas replace this part)
theta_clayton <- iTau(claytonCopula(), tau_target)
clayton_copula <- claytonCopula(param = theta_clayton, dim = 2)
df <- as.data.frame(rCopula(sample_size, clayton_copula))
colnames(df) <- paste0("u", 1:2)

# Generate survival times
df <- df %>%
  mutate(
    t1 = FAdist::qweibull3(u1, shape = nu1, scale = lambda1, thres = theta1),
    t2 = FAdist::qweibull3(u2, shape = nu2, scale = lambda2, thres = theta2)
    ) %>%
  mutate(
    x1 = t1,
    x2 = t2,
    delta1 = 1,  # without censoring
    delta2 = 1
  )

# Perform homogeneity tests (survival functions are the same)
akritas_related_samples(df)
sign_related_samples(df)
wilcoxon_signed_rank_related_samples(df)
weighted_logrank_related_samples(df)
ppw_related_samples(df)


# Censored data, two-sample case ----
sample_size <- 500 # sample size
tau_target <- 2/3 # Kendall tau
lambda1 <- 1 # scale parameter for t1
lambda2 <-1 # scale parameter for t2
nu1 <- 1 # shape parameter for t1
nu2 <- 1 # shape parameter for t2
theta1 <- 0 # location parameter for t1
theta2 <- 0 # location parameter for t2 
censorship_rate <- 0.15 # censorship rate

# Generate data from the copula (for other copulas replace this part)
theta_clayton <- iTau(claytonCopula(), tau_target)
clayton_copula <- claytonCopula(param = theta_clayton, dim = 2)
df <- as.data.frame(rCopula(sample_size, clayton_copula))
colnames(df) <- paste0("u", 1:2)

# Generate survival times
df <- df %>%
  mutate(
    t1 = FAdist::qweibull3(u1, shape = nu1, scale = lambda1, thres = theta1),
    t2 = FAdist::qweibull3(u2, shape = nu2, scale = lambda2, thres = theta2)
  ) %>%
  mutate(
    x1 = t1,
    x2 = t2,
    delta1 = 1, 
    delta2 = 1
  )

# Calculate parameter for specific censorship rate
lambda_min <- 0.000001
lambda_max <- 1000000

f <- function(lambda){
  mean(punif(df$t1, min = 0, max = lambda)) - censorship_rate
}
# Find solution f (lambda) = 0
lambda_cens <- uniroot(f, interval = c (lambda_min , lambda_max))$root 
lambda_cens

# Censoring 
c<- runif(sample_size, 0, lambda_cens)
df <- df %>%
  mutate(
    c = c,
    x1 = ifelse(t1 > c, c, t1),
    x2 = ifelse(t2 > c, c, t2),
    delta1 = ifelse(t1 > c, 0, 1),
    delta2 = ifelse(t2 > c, 0, 1))

# Perform homogeneity tests (survival functions are the same)
akritas_related_samples(df)
weighted_logrank_related_samples(df)
ppw_related_samples(df)

# Complete data, three-sample case ----
sample_size <- 500 # sample sizes
tau_target <- 2/3 # Kendall tau
lambda1 <- 1 # scale parameter for t1
lambda2 <- 1 # scale parameter for t2
lambda3 <- 1 # scale parameter for t3
nu1 <- 1 # shape parameter for t1
nu2 <- 1 # shape parameter for t2
nu3 <- 1 # shape parameter for t3
theta1 <- 0 # location parameter for t1
theta2 <- 0 # location parameter for t2 
theta3 <- 0 # location parameter for t3 

# Generate data from the copula
theta_clayton <- iTau(claytonCopula(), tau_target)
clayton_copula <- claytonCopula(param = theta_clayton, dim = 3)
df <- as.data.frame(rCopula(sample_size, clayton_copula))
colnames(df) <- paste0("u", 1:3)

# Generate survival times (Weibull distribution)
df <- df %>%
  mutate(
    t1 = FAdist::qweibull3(u1, shape = nu1, scale = lambda1, thres = theta1),
    t2 = FAdist::qweibull3(u2, shape = nu2, scale = lambda2, thres = theta2),
    t3 = FAdist::qweibull3(u3, shape = nu3, scale = lambda3, thres = theta3)
  ) %>%
  mutate(
    x1 = t1,
    x2 = t2,
    x3 = t3,
    delta1 = 1, # without censoring
    delta2 = 1,
    delta3 = 1
    )

# Perform homogeneity tests (survival functions are the same)
friedman_related_samples(df)
weighted_logrank_related_samples_more_than_two(df)


# Censored data, three-sample case ----
sample_size <- 500 # sample sizes
tau_target <- 2/3 # Kendall tau
lambda1 <- 1 # scale parameter for t1
lambda2 <- 1 # scale parameter for t2
lambda3 <- 1 # scale parameter for t3
nu1 <- 1 # shape parameter for t1
nu2 <- 1 # shape parameter for t2
nu3 <- 1 # shape parameter for t3
theta1 <- 0 # location parameter for t1
theta2 <- 0 # location parameter for t2 
theta3 <- 0 # location parameter for t3 
censorship_rate <- 0.15 # censorship rate

# Generate data from the copula
theta_clayton <- iTau(claytonCopula(), tau_target)
clayton_copula <- claytonCopula(param = theta_clayton, dim = 3)
df <- as.data.frame(rCopula(sample_size, clayton_copula))
colnames(df) <- paste0("u", 1:3)

# Generate survival times (Weibull distribution)
df <- df %>%
  mutate(
    t1 = FAdist::qweibull3(u1, shape = nu1, scale = lambda1, thres = theta1),
    t2 = FAdist::qweibull3(u2, shape = nu2, scale = lambda2, thres = theta2),
    t3 = FAdist::qweibull3(u3, shape = nu3, scale = lambda3, thres = theta3)
  )

# Calculate parameter for specific censorship rate
lambda_min <- 0.000001
lambda_max <- 1000000

f <- function(lambda){
  mean(punif(df$t1, min = 0, max = lambda)) - censorship_rate
}
# Find solution f (lambda) = 0
lambda_cens <- uniroot(f, interval = c (lambda_min , lambda_max))$root 
lambda_cens

# Censoring 
c<- runif(sample_size, 0, lambda_cens)
df <- df %>%
  mutate(
    c = c,
    x1 = ifelse(t1 > c, c, t1),
    x2 = ifelse(t2 > c, c, t2),
    x3 = ifelse(t3 > c, c, t3),
    delta1 = ifelse(t1 > c, 0, 1),
    delta2 = ifelse(t2 > c, 0, 1),
    delta3 = ifelse(t3 > c, 0, 1)
  )

# Perform homogeneity test (survival functions are the same)
weighted_logrank_related_samples_more_than_two(df)


