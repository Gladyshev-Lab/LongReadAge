

# Load necessary libraries
library(testit)
library(tidyverse)
library(data.table)
library(pbapply)
library(grid)
library(glue)
library(ggtext)


#################################
#                               #
#           Functions           #
#                               #
#################################

### construct reference
# Constructs a reference dataset containing linear age-methylation relationships.
# Returns a dataframe with numerical parameters describing age-methylation 
# relationship for each CpG site derived from linear regression.
# 
# Args:
#   reference_beta_input: A data frame where columns 4 onwards contain beta values for CpG sites.
#   reference_ages_input: A numeric vector containing the ages corresponding to the beta values.
# 
# Returns:
#   A data frame with columns for chromosome, start, end, Pearson correlation with age, 
#   range of methylation values, intercept, and coefficient of the linear model.

construct_reference <- function(reference_beta_input, reference_ages_input) {
  
  # Helper function to find correlation and linear model coefficients for a given CpG site
  find_info <- function(x) {
    
    pearson_r <- unname(cor.test(x, reference_ages_input)$estimate)
    range <- max(x) - min(x)
    
    # create df with age + methylation
    age_meth_df <- data.frame(
      age = reference_ages_input,
      meth = x
    )
    
    meth_predictor <- lm(meth ~ age, data = age_meth_df)
    
    intercept <- unname(meth_predictor$coefficients[1])
    coefficient <- unname(meth_predictor$coefficients[2])
    
    return(c(pearson_r, range, intercept, coefficient))
  }
  
  # Apply the helper function to each CpG site and construct the reference info data frame
  reference_info <- pbapply(reference_beta_input[, 4:ncol(reference_beta_input)], 1, find_info)
  reference_info <- as.data.frame(t(reference_info))
  
  colnames(reference_info) <- c("pearson_r", "range", "intercept", "coef")
  reference_info <- cbind(reference_beta_input[, 1:3], reference_info)
  colnames(reference_info)[1:3] <- c("chr", "start", "end")
  
  return(reference_info)
  
}


### generate age predictions
# Generates age predictions based on a reference dataset and sample methylation data.
# 
# Args:
#   reference_info_input: A data frame containing reference information about CpG sites.
#   sample_meth_vector: A data frame containing methylation data for the sample.
#   percentile_cutoff: A numeric value indicating the percentile cutoff for selecting CpG sites.
#   number_cutoff: A numeric value indicating the number cutoff for selecting CpG sites.
#   cutoff_mode: A string indicating the mode of cutoff ("percentile" or "number").
#   start: A numeric value indicating the starting age for prediction.
#   stop: A numeric value indicating the stopping age for prediction.
#   step: A numeric value indicating the step size for the age range.
# 
# Returns:
#   A numeric vector containing the predicted age, mean Pearson R for predictors, and number of intersecting CpG sites.

generate_age_predictions <- function(reference_info_input, sample_meth_vector, percentile_cutoff, number_cutoff, cutoff_mode, start, stop, step) {
  
  print(paste("Number of input CpGs", nrow(sample_meth_vector)))
  
  # intersect reference and sample methylation
  sample_meth_intersected <- inner_join(reference_info_input, sample_meth_vector, by = c("chr", "start", "end")) 
  print(paste("Number of intersected CpGs", nrow(sample_meth_intersected)))
  
  # sort by absolute value of pearson correlation of methylation with age in reference
  sample_meth_intersected <- sample_meth_intersected[order(abs(sample_meth_intersected$pearson_r), decreasing = T), ]
  
  if (cutoff_mode == "percentile") {
    
    # select top percentile_cutoff % of cpgs
    sample_meth_intersected <- sample_meth_intersected[1:(as.integer(nrow(sample_meth_intersected) * percentile_cutoff)), ]
    
  } else if (cutoff_mode == "number") {
  
    # select top # cutoff of cpgs
    sample_meth_intersected <- sample_meth_intersected[1:number_cutoff, ]
    
  }
  
  num_cpgs_intersected <- nrow(sample_meth_intersected)
  
  print(paste("Number of CpGs in model:", num_cpgs_intersected))
  
  if (num_cpgs_intersected == 0) {
    print("No intersecting CpGs")
  }
  
  # iterate through range of ages, see which one best predicts observed methylation pattern
  age_ranges <- seq(start, stop, step)
  prob_age <- numeric(length(age_ranges))
  
  observed_meth <- sample_meth_intersected[, 8]
  
  for (j in 1:length(age_ranges)) {
    
    predicted_meth <- numeric(num_cpgs_intersected)
    
    # compute predicted methylation given age
    for (m in 1:num_cpgs_intersected) {
      predicted_meth[m] <- age_ranges[j] * sample_meth_intersected$coef[m] + sample_meth_intersected$intercept[m]
    }
    
    # convert any predicted meth > 1 to 0.999, any predicted meth < 0 to 0.001
    predicted_meth[predicted_meth > 1] <- 0.999
    predicted_meth[predicted_meth < 0] <- 0.001
    
    # compute probability methylation associated with that age (sum log to avoid underflow errors)
    prob_age[j] <- sum(log(predicted_meth^observed_meth * (1 - predicted_meth)^(1 - observed_meth))) 
  }

  predicted_age <- age_ranges[which(prob_age == max(prob_age))]
  print(paste("Predicted age:", predicted_age))
  
  mean_r_predictors <- mean(abs(sample_meth_intersected$pearson_r))
  print(paste("Mean R predictors:", mean_r_predictors))
  
  # return: predicted age, mean pearson R for predictors, number of intersecting cpgs
  return(c(predicted_age, mean_r_predictors, num_cpgs_intersected))
}


