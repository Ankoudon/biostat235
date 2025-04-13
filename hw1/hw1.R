############### data ###############
control <- data.frame(
  S=c('F','F','F','M','M','M'),
  X=c(288,302,358,315,273,257),
  Yobs=c(293,306,367,320,279,254))

treatment <- data.frame(
  S=c('F','F','F','M','M','M'),
  X=c(285,299,366,334,266,251),
  Yobs=c(282,256,307,284,233,237))

data <- rbind(
  transform(control, Group='Control'),
  transform(treatment, Group='Treatment'))

############### 2.1 (a) ###############

# Wilcoxon Rank Sum test
wilcox.test(X ~ Group, data=data)
# Fisher's Exact test
table_sex <- table(data$S, data$Group)
fisher.test(table_sex)

############### 2.1 (b) ###############

# function (mean)

perm_test_exact <- function(values, st = "mean") {
  
  # Number of observations
  n <- length(values)
  # Check if n is even
  if (n %% 2 != 0) stop("n must be even.")
  
  # Group size
  group_size <- n / 2
  
  # Combinations
  combos <- combn(n, group_size)
  num_combos <- ncol(combos)
  
  # Observed difference
  original_labels <- 1:group_size
  
  if (st == "mean") {
    
    observed_diff <- mean(values[original_labels]) - mean(values[-original_labels])
    
    # Permutations
    perm_diffs <- apply(combos, 2, function(treatment_idx) {
      control_idx <- setdiff(1:n, treatment_idx)
      mean(values[treatment_idx]) - mean(values[control_idx])
    })
    
  } else if (st == "median") {
    
    observed_diff <- median(values[original_labels]) - median(values[-original_labels])
    
    # Permutations using median
    perm_diffs <- apply(combos, 2, function(treatment_idx) {
      control_idx <- setdiff(1:n, treatment_idx)
      median(values[treatment_idx]) - median(values[control_idx])
    })
    
  } else if (st == "var") {
    
    observed_diff <- var(values[original_labels]) - var(values[-original_labels])
    
    # Permutations using variance
    perm_diffs <- apply(combos, 2, function(treatment_idx) {
      control_idx <- setdiff(1:n, treatment_idx)
      var(values[treatment_idx]) - var(values[control_idx])
    })
    
  } else {
    
    stop("Invalid statistic type. Choose 'mean', 'median', or 'var'.")
    
  }
  
  # Calculate p-value
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  list(
    observed_difference = observed_diff,
    p_value = p_value)
}

# Observation
Y <- c(data$Yobs[data$Group == 'Control'], data$Yobs[data$Group == 'Treatment'])

perm_test_exact(Y, "mean")
perm_test_exact(Y, "median")
perm_test_exact(Y, "var")

############### 2.1 (c) ###############

# function (mean) with blocking
perm_test_blocked <- function(values, block, st = "mean") {
  
  if (length(values) != length(block)) stop("values and block must have same length.")
  n <- length(values)
  block <- as.factor(block)
  
  block_levels <- levels(block)
  
  block_combos <- lapply(block_levels, function(bl) {
    idx <- which(block == bl)
    k <- length(idx)
    if (k %% 2 != 0) stop(paste("Block", bl, "must have even number of observations."))
    combn(idx, k / 2, simplify = FALSE)
  })
  
  if (st == "mean") {
    
    observed_diff <- 0
    for (bl in block_levels) {
      idx <- which(block == bl)
      bl_values <- values[idx]
      group_size <- length(idx) / 2
      observed_diff <- observed_diff +
        mean(bl_values[1:group_size]) - mean(bl_values[-(1:group_size)])
    }
    
    combos_grid <- expand.grid(block_combos)
    
    perm_diffs <- apply(combos_grid, 1, function(row) {
      total_diff <- 0
      treatment_all <- numeric()
      control_all <- numeric()
      for (i in seq_along(block_levels)) {
        treatment_idx <- unlist(row[[i]])
        all_idx <- which(block == block_levels[i])
        control_idx <- setdiff(all_idx, treatment_idx)
        total_diff <- total_diff +
          mean(values[treatment_idx]) - mean(values[control_idx])
      }
      total_diff
    })
    
  } else if (st == "median") {
    
    observed_diff <- 0
    for (bl in block_levels) {
      idx <- which(block == bl)
      bl_values <- values[idx]
      group_size <- length(idx) / 2
      observed_diff <- observed_diff +
        median(bl_values[1:group_size]) - median(bl_values[-(1:group_size)])
    }
    
    combos_grid <- expand.grid(block_combos)
    
    perm_diffs <- apply(combos_grid, 1, function(row) {
      total_diff <- 0
      treatment_all <- numeric()
      control_all <- numeric()
      for (i in seq_along(block_levels)) {
        treatment_idx <- unlist(row[[i]])
        all_idx <- which(block == block_levels[i])
        control_idx <- setdiff(all_idx, treatment_idx)
        total_diff <- total_diff +
          median(values[treatment_idx]) - median(values[control_idx])
      }
      total_diff
    })
    
  } else if (st == "var") {
    
    observed_diff <- 0
    for (bl in block_levels) {
      idx <- which(block == bl)
      bl_values <- values[idx]
      group_size <- length(idx) / 2
      observed_diff <- observed_diff +
        var(bl_values[1:group_size]) - var(bl_values[-(1:group_size)])
    }
    
    combos_grid <- expand.grid(block_combos)
    
    perm_diffs <- apply(combos_grid, 1, function(row) {
      total_diff <- 0
      treatment_all <- numeric()
      control_all <- numeric()
      for (i in seq_along(block_levels)) {
        treatment_idx <- unlist(row[[i]])
        all_idx <- which(block == block_levels[i])
        control_idx <- setdiff(all_idx, treatment_idx)
        total_diff <- total_diff +
          var(values[treatment_idx]) - var(values[control_idx])
      }
      total_diff
    })
    
  } else {
    
    stop("Invalid statistic type. Choose 'mean', 'median', or 'var'.")
    
  }
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  list(
    observed_difference = observed_diff,
    p_value = p_value
  )
}


Y <-  c(data[data$Group == 'Control' & data$S == 'M', 'Yobs'],
        data[data$Group == 'Treatment'& data$S == 'M', 'Yobs'],
        data[data$Group == 'Control' & data$S == 'F', 'Yobs'],
        data[data$Group == 'Treatment'& data$S == 'F', 'Yobs'])

block <- c(rep('M', 3), rep('M', 3), rep('F', 3), rep('F', 3))

perm_test_blocked(Y, block, "mean")
perm_test_blocked(Y, block, "median")
perm_test_blocked(Y, block, "var")

############### 2.1 (d-b) ###############

Y <- c(data$Yobs[data$Group == 'Control'] - data$X[data$Group == 'Control'],
       data$Yobs[data$Group == 'Treatment'] - data$X[data$Group == 'Treatment'])

perm_test_exact(Y, "mean")
perm_test_exact(Y, "median")
perm_test_exact(Y, "var")

############### 2.1 (d-c) ###############

Y <-  c(
  data[data$Group == 'Control' & data$S == 'M', 'Yobs'] -
    data[data$Group == 'Control' & data$S == 'M', 'X'],
  data[data$Group == 'Treatment'& data$S == 'M', 'Yobs'] -
    data[data$Group == 'Treatment'& data$S == 'M', 'X'],
  data[data$Group == 'Control' & data$S == 'F', 'Yobs'] -
    data[data$Group == 'Control' & data$S == 'F', 'X'],
  data[data$Group == 'Treatment'& data$S == 'F', 'Yobs'] - 
    data[data$Group == 'Treatment'& data$S == 'F', 'X'])

block <- c(rep('M', 3), rep('M', 3), rep('F', 3), rep('F', 3))

perm_test_blocked(Y, block, "mean")
perm_test_blocked(Y, block, "median")
perm_test_blocked(Y, block, "var")

############### 2.2 (a) ###############

mu <- mean(treatment$Yobs) -  mean(control$Yobs) 
mu

############### 2.2 (b) ###############

var_treatment <- sd(treatment$Yobs)^2
var_treatment
var_control <- sd(control$Yobs)^2
var_control
var_fp <- var_treatment/6 + var_control/6
var_fp

ci_09 <- 1.645 * var_fp^(0.5)
c(mu - ci_09, mu + ci_09)

############### 2.2 (c) ###############

# Regression with W 
model2 <- lm(Yobs ~ Group, data = data)
summary(model2)


# Regression with W and X
model3 <- lm(Yobs ~ Group + X, data = data)
summary(model3)

############### 2.2 (d) ###############

gain_t <- treatment$Yobs - treatment$X
gain_c <- control$Yobs - control$X
mu <- mean(gain_t) -  mean(gain_c) 
mu

var_treatment <- sd(gain_t)^2
var_treatment
var_control <- sd(gain_c)^2
var_control
var_fp <- var_treatment/6 + var_control/6
var_fp

ci_09 <- 1.645 * var_fp^(0.5)
c(mu - ci_09, mu + ci_09)

