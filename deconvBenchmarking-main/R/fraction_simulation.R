##################### fraction simulation ##################
#' Simulate realistic cell-type fractions from beta distribution
#' @description Simulate cell-type fractions by utilizing the Cell-Type Proportion Distribution derived from scRNA.
#'    This approach is recommended when the scMeta contains a sufficient number of distinct samples to accurately capture the distribution
#'    and generate realistic proportions (n >= 10 is recommended).
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param n number of simulated bulk samples
#' @param colnames_of_sample column name that corresponds to cellType in scMeta
#' @param colnames_of_cellType column name that corresponds to sampleID in scMeta
#' @param fixed_cell_type a character denotes the target cell type for which we strive to faithfully preserve its distribution.
#'    It is recommended to set this parameter to the name of the malignant cell types. If left undefined, the function will automatically
#'    select the most abundant cell type as 'fixed_cell_type'.
#' @param min.frac minimum fraction in the simulated fraction, values below this threshold will be set to zero. Default = 0.01
#' @param showFractionPlot a logical variable determining whether to display simulated fraction distribution for the fixed_cell_type
#'
#' @return a matrix of simulated fraction, with samples in rows and cell_types in columns
#' @export
#'
#' @examples
#' \dontrun{
#' fracSimulator_Beta(scMeta, n = 100, colnames_of_sample = 'sampleID',
#'                    colnames_of_cellType = 'cell_type', fixed_cell_type = 'malignant')
#' }
fracSimulator_Beta<-function(scMeta,n,
                             colnames_of_sample = NA,
                             colnames_of_cellType = NA,
                             fixed_cell_type = NA,
                             min.frac = 0.01,
                             showFractionPlot = T){

  if(is.na(colnames_of_sample)|is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to sampleID and cell_type in scMeta')
  }

  ct_table = scMeta %>%
    dplyr::group_by(!!rlang::sym(colnames_of_sample), !!rlang::sym(colnames_of_cellType)) %>% # MODIFIED LINE
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(freq = n / sum(n))

  if(is.na(fixed_cell_type)){
    # Use the variable holding the actual column name for cell types
    cell_type_column_values <- scMeta[[colnames_of_cellType]] 
    fixed_cell_type <- names(table(cell_type_column_values))[which.max(table(cell_type_column_values))]
    message(paste0("No fixed_cell_type provided, automatically selected: ", fixed_cell_type))
  }
  
  message(paste0("Processing fixed_cell_type: ", fixed_cell_type))
  
  # 1. Get the original observed frequencies for the fixed cell type
  observed_freq_for_fixed <- ct_table[['freq']][ct_table[[colnames_of_cellType]] == fixed_cell_type]
  
  # Use x_fixed for modifications needed for fitting beta distribution
  x_fixed <- observed_freq_for_fixed 
  
  # Debugging prints (optional)
  # print(paste0("Number of observations for ", fixed_cell_type, " (variable x_fixed initial): ", length(x_fixed)))
  # if(length(x_fixed) > 0) print(summary(x_fixed))
  
  # 2. Logic to determine fixed_frac based on x_fixed (derived from observed_freq_for_fixed)
  if(length(x_fixed) > 1) {
    # Adjust values to be strictly within (0, 1) for beta distribution
    epsilon <- 1e-6 # A small value to nudge away from 0 and 1
    x_fixed_for_fitting <- x_fixed # Create a copy for fitting modifications
    x_fixed_for_fitting[x_fixed_for_fitting >= 1] <- 1 - epsilon
    x_fixed_for_fitting[x_fixed_for_fitting <= 0] <- epsilon
    # Ensure all are strictly within (0,1) after nudging
    x_fixed_for_fitting <- pmax(epsilon, pmin(x_fixed_for_fitting, 1 - epsilon)) 
    
    fit_results_fixed <- tryCatch({
      fitdistrplus::fitdist(x_fixed_for_fitting, "beta", lower = c(0, 0)) 
    }, error = function(e) {
      warning(paste("Could not fit beta distribution for fixed_cell_type '", fixed_cell_type, 
                    "'. Error: ", e$message, 
                    ". Attempted fit on data: ", paste(round(x_fixed_for_fitting,3), collapse=", ")), call. = FALSE)
      return(NULL)
    })
    
    if (!is.null(fit_results_fixed) && !is.null(fit_results_fixed$estimate) && all(!is.na(fit_results_fixed$estimate))) {
      fixed_frac <- rbeta(n, fit_results_fixed$estimate[1], fit_results_fixed$estimate[2])
    } else {
      warning(paste("Using resampling from original observed frequencies for fixed_cell_type '", fixed_cell_type, 
                    "' due to fitting issues (e.g., too few unique values, fit failure, or NA estimates)."), call. = FALSE)
      # Resample from the original observed_freq_for_fixed, not the nudged version
      fixed_frac <- sample(observed_freq_for_fixed, size = n, replace = TRUE)
    }
  } else { 
    # This block handles length(x_fixed) <= 1 (i.e. 0 or 1 elements)
    warning(paste("Not enough data points (", length(x_fixed), ") to fit beta distribution for fixed_cell_type '", 
                  fixed_cell_type, "'. Using observed value(s) directly or a default."), call. = FALSE)
    if (length(x_fixed) == 1) { # If exactly one value, replicate it
      fixed_frac <- rep(x_fixed[1], n)
    } else { # length is 0 (this was the case for 'Mal' as per your warnings)
      warning(paste0("Fixed cell type '", fixed_cell_type, "' had no valid frequency data in ct_table for fitting. Assigning a small default fraction (0.01)."), call. = FALSE)
      fixed_frac <- rep(0.01, n) 
    }
  }
  # At this point, 'fixed_frac' is defined.
  # And 'observed_freq_for_fixed' holds the original unmodified frequencies.
  
  # 3. Define and call quick_hist
  quick_hist <- function(observed_data_for_plot, simulated_data_for_plot, ct_name_for_plot) {
    par(mfrow=c(1,2))
    if (length(observed_data_for_plot) > 0) {
      hist(observed_data_for_plot, 
           main = paste('Original observed frac for\n', ct_name_for_plot), 
           xlab="Frequency", breaks=seq(0,1,by=0.1), xlim=c(0,1))
    } else {
      plot(0, type="n", xlim=c(0,1), ylim=c(0,1), 
           xlab="Frequency", ylab="Count", 
           main = paste('Original observed frac for\n', ct_name_for_plot))
      text(0.5, 0.5, paste("No observed data for\n", ct_name_for_plot))
    }
    hist(simulated_data_for_plot, 
         main = paste('Simulated frac for\n', ct_name_for_plot), 
         xlab="Frequency", breaks=seq(0,1,by=0.1), xlim=c(0,1))
  }
  
  if(showFractionPlot==T){ # showFractionPlot defaults to T
    quick_hist(observed_data_for_plot = observed_freq_for_fixed,  # Use the original, unmodified observations
               simulated_data_for_plot = fixed_frac, 
               ct_name_for_plot = fixed_cell_type)
    par(mfrow=c(1,1)) # Reset plotting layout
  }

  # for other cell types, simulate cell fraction based on beta distribution, and scaled to the remaining frac
  # Modify quick_rbeta function similarly
  quick_rbeta <- function(current_cell_type_for_loop) {
    x_other <- ct_table[['freq']][ct_table[[colnames_of_cellType]] == current_cell_type_for_loop]
    
    # Debugging prints (optional)
    print(paste0("Processing other_cell_type: ", current_cell_type_for_loop, ", N_obs: ", length(x_other)))
    if(length(x_other) > 0) print(summary(x_other))
    
    if(length(x_other) > 1) {
      epsilon <- 1e-6
      x_other[x_other >= 1] <- 1 - epsilon
      x_other[x_other <= 0] <- epsilon
      x_other <- pmax(pmin(x_other, 1 - epsilon), epsilon)
      
      fit_x_other <- tryCatch({
        fitdistrplus::fitdist(x_other, "beta", lower = c(0, 0))
      }, error = function(e) {
        warning(paste("Could not fit beta for other_cell_type '", current_cell_type_for_loop, 
                      "'. Error: ", e$message,
                      ". Observed frequencies: ", paste(round(x_other,3), collapse=", ")), call. = FALSE)
        return(NULL)
      })
      
      if(!is.null(fit_x_other) && !is.null(fit_x_other$estimate) && all(!is.na(fit_x_other$estimate))){
        return(rbeta(n, fit_x_other$estimate[1], fit_x_other$estimate[2]))
      } else {
        warning(paste("Using resampling from observed frequencies for other_cell_type '", current_cell_type_for_loop, 
                      "' due to fitting issues."), call. = FALSE)
        if (length(x_other) > 0) {
          original_x_for_other <- ct_table[['freq']][ct_table[[colnames_of_cellType]] == current_cell_type_for_loop]
          return(sample(original_x_for_other, size = n, replace = TRUE))
        }
        return(rep(0, n)) # Fallback if no data and no fit
      }
    } else {
      warning(paste("Not enough data (", length(x_other), ") to fit beta for other_cell_type '", 
                    current_cell_type_for_loop, "'. Using observed frequencies or 0."), call. = FALSE)
      if (length(x_other) == 1) return(rep(x_other[1], n))
      return(rep(0, n)) # Fallback for no data
    }
  }

  all_cell_types_in_input_scMeta <- unique(scMeta[[colnames_of_cellType]]) # Reminder: colnames_of_cellType is "cell_type"
  other_cell_types <- all_cell_types_in_input_scMeta[all_cell_types_in_input_scMeta != fixed_cell_type]
  
  # Make sure the quick_rbeta function is defined here with all the robust fitting logic
  # (length checks, data nudging, tryCatch, fallbacks similar to fixed_cell_type)
  # It should consistently return a vector of length 'n', e.g., rep(0, n) on failure.
  quick_rbeta <- function(current_cell_type_for_loop) {
    x_other <- ct_table[['freq']][ct_table[[colnames_of_cellType]] == current_cell_type_for_loop]
    
    if(length(x_other) > 1) {
      epsilon <- 1e-6
      data_for_fitting_other <- x_other # Create a copy for modification
      data_for_fitting_other[data_for_fitting_other >= 1] <- 1 - epsilon
      data_for_fitting_other[data_for_fitting_other <= 0] <- epsilon
      data_for_fitting_other <- pmax(epsilon, pmin(data_for_fitting_other, 1 - epsilon))
      
      fit_results_other <- tryCatch({
        fitdistrplus::fitdist(data_for_fitting_other, "beta", lower = c(0, 0))
      }, error = function(e) {
        warning(paste("Could not fit beta for other_cell_type '", current_cell_type_for_loop, 
                      "'. Error: ", e$message,
                      ". Observed frequencies: ", paste(round(data_for_fitting_other,3), collapse=", ")), call. = FALSE)
        return(NULL)
      })
      
      if(!is.null(fit_results_other) && !is.null(fit_results_other$estimate) && all(!is.na(fit_results_other$estimate))){
        return(rbeta(n, fit_results_other$estimate[1], fit_results_other$estimate[2]))
      } else {
        warning(paste("Using resampling from observed for other_cell_type '", current_cell_type_for_loop, 
                      "' due to fitting issues."), call. = FALSE)
        if (length(x_other) > 0) { # Use original x_other for resampling
          return(sample(x_other, size = n, replace = TRUE))
        }
        return(rep(0, n)) # Fallback if no data and no fit
      }
    } else {
      warning(paste("Not enough data (", length(x_other), ") to fit beta for other_cell_type '", 
                    current_cell_type_for_loop, "'. Using observed or 0."), call. = FALSE)
      if (length(x_other) == 1) return(rep(x_other[1], n))
      return(rep(0, n)) # Fallback for no data
    }
  }
  # End of quick_rbeta definition
  
  cat("\n--- DEBUG: INSIDE fracSimulator_Beta, before processing other_cell_types ---\n")
  cat("Value of fixed_cell_type:", fixed_cell_type, "\n")
  cat("Unique cell types in scMeta input to fracSimulator_Beta (all_cell_types_in_input_scMeta):\n")
  print(all_cell_types_in_input_scMeta)
  cat("Calculated other_cell_types:\n")
  print(other_cell_types)
  cat("Length of other_cell_types:", length(other_cell_types), "\n")
  
  if (length(other_cell_types) > 0) {
    cat("DEBUG fracSimulator_Beta: Entering block for length(other_cell_types) > 0\n") # Add this
    
    flexible_frac_list <- lapply(other_cell_types, quick_rbeta)
    flexible_frac <- do.call(cbind, flexible_frac_list)
    
    # Ensure flexible_frac is a matrix, even if only one other_cell_type
    if (!is.matrix(flexible_frac) && length(flexible_frac) == n && length(other_cell_types) == 1) {
      flexible_frac <- matrix(flexible_frac, ncol = 1)
    }
    
    if (!is.null(flexible_frac) && is.matrix(flexible_frac) && ncol(flexible_frac) > 0) {
      remaining_frac <- 1 - fixed_frac 
      # Ensure remaining_frac is not negative (can happen if fixed_frac due to rbeta is > 1, though unlikely with fitted params)
      remaining_frac[remaining_frac < 0] <- 0
      
      sum_flexible_frac <- rowSums(flexible_frac)
      scaling_factor <- rep(0, n) # Initialize
      
      # Calculate scaling factor only for rows where sum_flexible_frac is positive
      # and remaining_frac is also positive (to avoid 0/0 or x/0 if remaining_frac is also 0)
      can_scale_indices <- sum_flexible_frac > 1e-9 & remaining_frac > 1e-9
      if (any(can_scale_indices)) {
        scaling_factor[can_scale_indices] <- remaining_frac[can_scale_indices] / sum_flexible_frac[can_scale_indices]
      }
      
      flexible_frac <- sweep(flexible_frac, 1, scaling_factor, '*')
      
      # For rows where sum_flexible_frac was 0, but remaining_frac > 0,
      # distribute remaining_frac equally among other_cell_types for those rows.
      distribute_equally_indices <- (sum_flexible_frac <= 1e-9) & (remaining_frac > 1e-9)
      if (any(distribute_equally_indices) && length(other_cell_types) > 0) {
        warning("Some rows had zero sum for flexible_frac but non-zero remaining_frac. Distributing equally.", call. = FALSE)
        num_other_types_for_distrib <- max(1, length(other_cell_types))
        for(idx in which(distribute_equally_indices)){
          flexible_frac[idx, ] <- remaining_frac[idx] / num_other_types_for_distrib
        }
      }
      
      simulated_frac <- cbind(fixed_frac, flexible_frac)
      colnames(simulated_frac) <- c(fixed_cell_type, other_cell_types)
    } else {
      # This case means other_cell_types existed, but all quick_rbeta calls might have failed to produce usable vectors,
      # or cbind resulted in something not a matrix.
      warning("Processing of other_cell_types did not yield valid matrix for flexible_frac. Using only fixed_cell_type.", call. = FALSE)
      simulated_frac <- matrix(1.0, nrow = n, ncol = 1) # Make fixed_cell_type 100%
      colnames(simulated_frac) <- fixed_cell_type
    }
  } else {
    # No other_cell_types, fixed_cell_type is the only one.
    cat("DEBUG fracSimulator_Beta: Entering block for length(other_cell_types) == 0\n") # Add this
    message(paste0("Only fixed_cell_type '", fixed_cell_type, "' is present/selected. Its fraction will be 100%."))
    simulated_frac <- matrix(1.0, nrow = n, ncol = 1)
    colnames(simulated_frac) <- fixed_cell_type
  }
  
  # Final normalization and cleaning
  simulated_frac[is.na(simulated_frac)] <- 0 # Catch any NAs
  simulated_frac[simulated_frac < min.frac] <- 0 # Apply min.frac threshold
  
  # Normalize rows to sum to 1. Handle rows that sum to 0 after thresholding.
  row_sums_final <- rowSums(simulated_frac)
  for (i in seq_len(nrow(simulated_frac))) {
    if (row_sums_final[i] < 1e-9) { # If sum is effectively zero
      warning(paste0("Row ", i, " of simulated fractions summed to 0 after min.frac. Distributing equally among all ", ncol(simulated_frac), " considered types."), call. = FALSE)
      if (ncol(simulated_frac) > 0) {
        simulated_frac[i, ] <- 1 / ncol(simulated_frac)
      } else {
        # This case should ideally not be reached if fixed_cell_type always results in a column
        simulated_frac[i, ] <- 0 # Or handle as error
      }
    } else {
      simulated_frac[i, ] <- simulated_frac[i, ] / row_sums_final[i]
    }
  }
  
  return(simulated_frac)
}


#' Simulate cell-type fractions from dirichlet distribution
#'
#' @param frac_prior A named numeric vector representing the fraction prior. The simulated fraction will fluctuate around the values specified in this vector.
#' @param n number of random fractions to generate
#' @param dispersion_par a numeric value determine the dispersion level of the simulated fractions. With lower value indicating higher dispersion level.
#' @param min.frac minimum fraction in the simulated fraction, values below this threshold will be set to zero
#' @param showFractionPlot a logical variable determining whether to display the boxplot of simulated fractions. It is recommended to set this parameter to TRUE to
#'    visualize the effect dispersion_par and choose the optimal value of it
#'
#' @return a matrix of simulated fraction, with samples in rows and cell_types in columns
#' @export
#'
#' @examples
#' frac_prior = c(100,200,300,400)
#' names(frac_prior) = paste0('cellType',seq(1,4))
#' fracSimulator_Dirichlet(frac_prior, n = 10, dispersion_par = 0.01)

fracSimulator_Dirichlet<-function(frac_prior,n, dispersion_par = 0.1, min.frac = 0, showFractionPlot=T){

  simulated_frac = gtools::rdirichlet(n,frac_prior*dispersion_par)
  colnames(simulated_frac) = names(frac_prior)

  simulated_frac[simulated_frac<min.frac] = 0
  simulated_frac = sweep(simulated_frac,1,rowSums(simulated_frac),'/')

  if(showFractionPlot){
    boxplot(simulated_frac)
  }

  return(simulated_frac)
}


#' Simulate cell-type fractions using Generator() function from Favilaco et al.
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param nbulk number of simulated bulk samples
#' @param pool.size number of cells to aggregate for each simulated bulk sample
#' @param min.percentage minimum percentage of cellType fraction to generate in fraction simulation. Default = 1
#' @param max.percentage maximum percentage of cellType fraction to generate in fraction simulation. Default = 99
#' @param seed a single value. Default = 24
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
fracSimulator_favilaco <- function(scExpr, scMeta, colnames_of_cellType = NA, nbulk = 100, pool.size = 100,
                                   min.percentage = 1, max.percentage = 99, seed = 24){

  # this function incorporated the Generator() function from https://github.com/favilaco/deconv_benchmark/blob/master/helper_functions.R
  # the Generator() function expects the 'phenoData' argument to contain two columns named 'cellID' and 'cellType'
  # Furthermore, this function does not require a pre-defined proportion matrix, it will generate and provide the simulated fraction matrix autmatically

  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  phenoData <- scMeta %>% mutate_('cellType'=colnames_of_cellType)
  phenoData$cellID = rownames(phenoData)

  v = Generator(scExpr,phenoData,Num.mixtures = nbulk, pool.size , min.percentage , max.percentage , seed )
  simulated_frac = as.matrix(t(v$P))
  rownames(simulated_frac) = NULL

  return(simulated_frac)
}


#' Simulate cell-type fractions using generateBulk_norep() function from SCDC package
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param disease indicate the health condition of subjects
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.
#' @param nbulk number of simulated bulk samples
#' @param samplewithRep logical, randomly sample single cells with replacement. Default is T
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
fracSimulator_SCDC = function(scExpr,scMeta,colnames_of_cellType = NA, colnames_of_sample = NA, disease = NULL, ct.sub = NULL,
                              nbulk = 10, samplewithRep = T, prop_mat = NULL){
  require(SCDC,quietly = T)
  eset = Biobase::ExpressionSet(assayData = scExpr,phenoData = new("AnnotatedDataFrame", data = scMeta))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.na(colnames_of_sample)){
    stop('please provide column name that corresponds to sampleID in scMeta')
  }

  if(is.null(ct.sub)){
    ct.sub = unique(scMeta[,colnames_of_cellType])
  }
  v = SCDC::generateBulk_norep(eset,ct.varname=colnames_of_cellType,sample = colnames_of_sample,disease, ct.sub, prop_mat, nbulk, samplewithRep)

  frac = v$true_p
  rownames(frac) = NULL

  return(frac)
}

