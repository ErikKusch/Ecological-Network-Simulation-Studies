#' Network Inference Behaviour Validation
#'
#' Logistic Bayesian models of network inference approach behaviour. Runs models per link sign (positive, negative, absent) to identify the effects of underlying link magnitude and environmental niche preference differences.
#'
#' @param model_df A data frame tallying comparisons of inferred and underlying true networks. Rows index pairwise links. Consists of the following columns:
#'  - Magnitude ... Numeric. Absolute strength/magnitude of link in the true network.
#'  - SignTrue ... Numeric - -1, 0, or 1. Sign of link in true network. 0 denotes absence of a link.
#'  - SignInferred ... Numeric - -1, 0, or 1. Sign of link in inferred network. 0 denotes absence of a link.
#'  - EnvDiff ... Numeric. Difference in environmental niche preference of linked species.
#' See examples for a workflow of how to create this data frame.
#' @param mode Character. Either "Inference" or "Detection". Mode of behaviour validation:
#'  - Inference ... likelihood of assignment of a link irrespective of whether it is correctly assigned
#'  - Detection ... likelihood of correct inference of individual links
#' @param nWarmup Integer. Number of warmup iterations for Bayesian models.
#' @param nSamples Integer. Number of total iterations for Bayesian models.
#' @param nChains Integer. Number of markov chains to run.
#' @param nCores Integer. Number of cores to utilise.
#' @param seed Numeric. Seed for random processes.
#'
#' @return A list with three elements:
#'  - Positive ... brmsfit object describing the Bayesian model run to identify either inference or detection behaviour for positive inferred links
#'  - Negative ... brmsfit object describing the Bayesian model run to identify either inference or detection behaviour for negative inferred links
#'  - Absent ... brmsfit object describing the Bayesian model run to identify either inference or detection behaviour for absent inferred links
#'
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#'
#' @importFrom brms brm
#' @importFrom brms bernoulli
#'
#' @examples
#' ## loading data ---
#' data("Network_mat")
#' data("Inferred_mat")
#' data("Niches_vec")
#'
#' matrices_ls <- list(True = Network_mat, Inferred = Inferred_mat)
#'
#' ### Calculate environmental differences
#' SPTrait_df <- data.frame(Niches_vec)
#' SPTrait_df$Species <- rownames(SPTrait_df)
#' colnames(SPTrait_df) <- c("Trait", "Species")
#' SPTrait_mat <- abs(outer(SPTrait_df$Trait, SPTrait_df$Trait, "-"))
#' colnames(SPTrait_mat) <- rownames(SPTrait_mat) <- colnames(SPTrait_df$Species)
#'
#' ## Collate data frame for behaviour analysis ---
#' model_df <- data.frame(
#'   Magnitude = abs(as.vector(matrices_ls[[1]])),
#'   SignTrue = as.vector(sign(matrices_ls[[1]])),
#'   SignInferred = as.vector(sign(matrices_ls[[2]])),
#'   EnvDiff = as.vector(SPTrait_mat)
#' )
#' model_df <- na.omit(model_df)
#' \dontrun{
#' (Detection <- Val.Behaviour(model_df, mode = "Detection"))
#' (Inference <- Val.Behaviour(model_df, mode = "Inference"))
#' }
#'
#' @export
Val.Behaviour <- function(model_df, mode = "Detection", nWarmup = 3e3, nSamples = 1e4, nChains = 4, nCores = 4, seed = 42) {
  model_df$CorrectPos <- model_df$SignInferred + model_df$SignTrue == 2
  model_df$CorrectNeg <- model_df$SignInferred + model_df$SignTrue == -2
  model_df$CorrectAbs <- (model_df$SignInferred == 0) + (model_df$SignTrue == 0) == 2

  if (mode == "Detection") {
    Bayes_Model_Positive <- brm(
      formula = CorrectPos ~ Magnitude * EnvDiff,
      data = model_df,
      family = bernoulli(link = "logit"),
      warmup = nWarmup,
      iter = nSamples,
      chains = nChains,
      cores = nCores,
      seed = seed
    )

    Bayes_Model_Negative <- brms::brm(
      formula = CorrectNeg ~ Magnitude * EnvDiff,
      data = model_df,
      family = bernoulli(link = "logit"),
      warmup = nWarmup,
      iter = nSamples,
      chains = nChains,
      cores = nCores,
      seed = seed
    )

    Bayes_Model_Absent <- brms::brm(
      formula = CorrectAbs ~ Magnitude * EnvDiff,
      data = model_df,
      family = bernoulli(link = "logit"),
      warmup = nWarmup,
      iter = nSamples,
      chains = nChains,
      cores = nCores,
      seed = seed
    )
  }
  if (mode == "Inference") {
    run_df <- model_df
    run_df$SignInferred[run_df$SignInferred != 1] <- 0
    Bayes_Model_Positive <- brm(
      formula = SignInferred ~ Magnitude * EnvDiff,
      data = run_df,
      family = bernoulli(link = "logit"),
      warmup = nWarmup,
      iter = nSamples,
      chains = nChains,
      cores = nCores,
      seed = seed
    )

    run_df <- model_df
    run_df$SignInferred[run_df$SignInferred != -1] <- 0
    run_df$SignInferred <- abs(run_df$SignInferred)
    Bayes_Model_Negative <- brm(
      formula = SignInferred ~ Magnitude * EnvDiff,
      data = run_df,
      family = bernoulli(link = "logit"),
      warmup = nWarmup,
      iter = nSamples,
      chains = nChains,
      cores = nCores,
      seed = seed
    )

    run_df <- model_df
    run_df$SignInferred[run_df$SignInferred != 0] <- 99
    run_df$SignInferred[run_df$SignInferred == 0] <- 1
    run_df$SignInferred[run_df$SignInferred == 99] <- 0
    Bayes_Model_Absent <- brm(
      formula = SignInferred ~ Magnitude * EnvDiff,
      data = run_df,
      family = bernoulli(link = "logit"),
      warmup = nWarmup,
      iter = nSamples,
      chains = nChains,
      cores = 1,
      seed = 42
    )
  }

  Return_ls <- list(
    Positive = Bayes_Model_Positive,
    Negative = Bayes_Model_Negative,
    Absent = Bayes_Model_Absent
  )
  return(Return_ls)
}
