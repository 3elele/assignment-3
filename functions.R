#----------2a
# Difference in the proportion of cases with a specific value between two groups.
#
# ARGUMENTS:
# d: a data frame or tibble
# var: the name of a column of d containing the dependent variable, provided as a string
# grouping_var: the name of a column of d containing a grouping variable, provided as a string
# group1: the value of grouping_var that corresponds to the first group
# group2: the value of grouping_var that corresponds to the second group
# target_value: the value of var that will be counted
#
# RETURN VALUE:
# The percentage of cases in which `var` was equal to `target_value` for the first group,
# minus the percentage of cases in which `var` was equal to `target_value` for the
# second group.
#
difference_in_proportion <- function(d, var, grouping_var, group1, group2,
                                     target_value) {
  percentageGroup1 <- nrow(dplyr::filter(d, UQ(as.name(grouping_var))==group1, UQ(as.name(var))==target_value)) / nrow(dplyr::filter(d, UQ(as.name(grouping_var))==group1))
  percentageGroup2 <- nrow(dplyr::filter(d, UQ(as.name(grouping_var))==group2, UQ(as.name(var))==target_value)) / nrow(dplyr::filter(d, UQ(as.name(grouping_var))==group2))
  return(percentageGroup1-percentageGroup2)
}

#----------2b
# Randomize the order of a column.
#
# ARGUMENTS:
# d: a data frame or tibble
# var: the name of a column of d containing the variable to randomize,
#      provided as a string
#
# RETURN VALUE:
# A data frame or tibble exactly the same as d, except with the order of
# var permuted randomly.
#
randomize <- function(d, var) {
  d[[var]] <- sample(d[[var]])
  return(d)
}
# Analize greater and smaller values of an observed one.
#
# ARGUMENTS:
# p: a data list
# 
# RETURN VALUE:
# the proportion of observations that are greater than or less than
# the observed value of the test statistic
# plus a slight adjustment to resolve the inaccuracy that would exist
#
permutation_pvalue_right <- function(p) {
  n_above <- sum(p$permuted >= p$observed)
  n_samples <- length(p$permuted)
  return((n_above + 1)/(n_samples + 1))
}
permutation_pvalue_left <- function(p) {
  n_below <- sum(p$permuted <= p$observed)
  n_samples <- length(p$permuted)
  return((n_below + 1)/(n_samples + 1))
}

# Perform a permutation test for two groups.
#
# ARGUMENTS:
# d: a data frame or tibble
# var: the name of the column in d on which the test statistic will be calculated,
#      provided as a string
# grouping_var: the name of the column in d which gives the grouping
# group1: the value of grouping_var corresponding to the first group
# group2: the value of grouping_var corresponding to the second group
# statistic: a function yielding a test statistic, which takes as input
#            a data frame, the name of a variable on which to calculate the
#            test statistic, the name of a grouping variable, the value of
#            the grouping variable corresponding to the first group, and
#            the value of the grouping variable corresponding to the second
#            group
# n_samples: the number of permutation samples to draw (default: 9999)
#
# RETURN VALUE:
#
# A list containing two elements:
#
#  - observed: the value of statistic() in d
#  - permuted: a vector containing the values of statistic() under n_samples
#              permutations
#
permutation_twogroups <- function(d, var, grouping_var, group1, group2, statistic,
                                  n_samples=9999, ...) {
  observed_statistic <- statistic(d, var, grouping_var, group1, group2, ...)
  permutation_statistics <- rep(0, n_samples)
  for (i in 1:n_samples) {
    permutations <- randomize(d, var)
    permutation_statistics[i] <- statistic(permutations, var, grouping_var, group1, group2, ...) 
  }
  result <- list(observed=observed_statistic,
                 permuted=permutation_statistics)
  return(result)
}

#----------3b
# Perform a permutation test for two implicit groups of binary trials summarized
# by the number of "successes" and trials for each of the two groups, using the
# difference in the proportion of successes (group 1 minus group 2) as a test
# statistic.
#
# ARGUMENTS:
# k1: the number of "successes" (i.e., observations of one of the two types) in group 1
# k2: the number of "successes" in group 2
# n1: the total number of trials in group 1
# n2: the total number of trials in group 2
# n_samples: the number of permutations (defaults to 9999)
#
# RETURN VALUE:
#
# A list containing two elements:
#
#  - observed: the value of statistic() in d
#  - permuted: a vector containing the values of statistic() under n_samples
#              permutations
#
permtest_difference_in_props <- function(k1, k2, n1, n2, n_samples=9999) {
  # Create a set of observations with exactly k1 and k2 successes
  obs_1 <- c(rep(TRUE, k1), rep(FALSE, n1 - k1)) # Individual observations from group 1
  obs_2 <- c(rep(TRUE, k2), rep(FALSE, n2 - k2)) # Individual observations from group 2
  observations <- c(obs_1, obs_2)
  # Permute this set of observations n_samples times, saving the result in a
  # matrix
  rep_observations <- matrix(rep(observations, n_samples), n1 + n2)
  perm_observations <- apply(rep_observations, 2, sample, n1 + n2)
  # Generate the proportions in the two groups amongst the permuted observations.
  # Tricks: mean() of a TRUE/FALSE variable is the proportion "TRUE";
  # instead of having explicit "Group" labels that we hold fixed, we just hold fixed
  # that the first n1 rows are "Group 1" and the remaining n2 rows are "Group 2",
  # which amounts to the same thing, and we generate the two percentages directly.
  props_group1 <- colMeans(perm_observations[1:n1,])
  props_group2 <- colMeans(perm_observations[(n1+1):(n1+n2),])
  test_stats <- props_group1 - props_group2
  return(list(observed=((k1/n1) - (k2/n2)), permuted=test_stats))
}

v_pdp_pvalue_right <- function(k1_vec, k2_vec, n1, n2, n_samples=9999) {
  result <- rep(NA, length(k1_vec))
  for (i in 1:length(k1_vec)) {
    # [YOUR CODE HERE: APPLY permtest_difference_in_props WITH THE i'TH VALUES
    #  OF k1_vec AND OF k2_vec AS THE FIRST TWO ARGUMENTS, AND STORE THE
    #  RESULT AS THE i'TH VALUE OF result]
    result[i] <- permtest_difference_in_props(k1_vec[i],k2_vec[i],n1,n2,n_samples) %>% permutation_pvalue_right
      
  }
  return(result)
}

stat_power <- function(d, n_samples, t) {
  pvals <- d %>% dplyr::select(Permutation_Pvalue)
  return(1 - sum(pvals > t)/n_samples)
}

#----------4
# Pearson's cumulative test statistic
#
# ARGUMENTS:
# d: a data frame or tibble
# var: the name of a column of d, provided as a string
# grouping_var: the name of another column of d, provided as a string
#
# RETURN VALUE:
# The value of Pearson's cumulative test statistic, calculated over all possible
# combinations of `var` and `grouping_var`
pearson_x2_stat <- function(d, var, grouping_var) {
  values <- unique(d[[var]])
  groups <- unique(d[[grouping_var]])
  result <- 0
  for (v in values) {
    prop_overall <- mean(d[[var]] == v)
    for (g in groups) {
      d_g <- dplyr::filter(d, get(grouping_var) == g)
      prop_g <- mean(d_g[[var]] == v)
      result <- result + ((prop_g - prop_overall)^2)*(nrow(d_g)/prop_overall)
    }
  }
  return(result)
}


# Perform a permutation test.
#
# ARGUMENTS:
# d: a data frame or tibble
# var_to_permute: the name of the column of d to permute
# statistic: a function yielding a test statistic, which takes a table as input
# n_samples: the number of permutation samples to draw (default: 9999)
# ... : further arguments that will be passed to `statistic`
#
# RETURN VALUE:
# 
# A list containing two elements:
#
#  - observed: the value of statistic() in d
#  - permuted: a vector containing the values of statistic() under n_samples
#              permutations
#
permutation_test <- function(d, var_to_permute, statistic, n_samples=9999, ...) {
  observed_statistic <- statistic(d, ...)
  permutation_statistics <- rep(0, n_samples)
  for (i in 1:n_samples) {
    d_perm <- randomize(d, var_to_permute)
    permutation_statistics[i] <- statistic(d_perm, ...)
  }
  result <- list(observed=observed_statistic,
                 permuted=permutation_statistics)
  return(result)
}
