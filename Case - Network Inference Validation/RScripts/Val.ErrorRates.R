#' Network Inference Error Rates Validation
#'
#' Calculate accuracy of network inference by sign of underlying true link.
#'
#' @param Network1 matrix object. First network in comparison.
#' @param Network2 matrix object. Second network in comparison.
#' @param NetworkType Character. Whether to create an interaction or association matrix. In an interaction matrix, columns act on rows. In an association matrix, entries are mirrored around the diagonal. Possible values: "Association" or "Interaction".
#'
#' @return A dataframe containing percentages and error rate types:
#'    - TP ... true positive associations; how many of the inferred positive associations are correctly identified as such? - NaN means no positives in the inferred matrix
#'    - TN ... true negative associations; how many of the inferred negative associations are correctly identified as such? - NaN means no negatives in the inferred matrix
#'    - TA ... true absent; how many of the inferred absent interactions are correctly identified as such? - NaN means no absents in the inferred matrix
#'    - FP ... falsely inferred positive; how many of the inferred positive associations are incorrectly identified as such?
#'    - FN ... falsely inferred negative; how many of the inferred negative associations are incorrectly identified as such?
#'    - FA ... false absent; how many of the inferred absent interactions are correctly identified as such?
#'    - MP ... true positive links which were not inferred; how many of the true positive associations were not inferred?
#'    - MN ... true negative links which were not inferred; how many of the true negative associations were not inferred?
#'    - MA ... true absent links which were not inferred; how many of the true absent associations were not inferred?
#'
#'
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#'
#' @examples
#' data("Network_mat")
#' data("Inferred_mat")
#' Val.ErrorRates(Network1 = Network_mat, Network2 = Inferred_mat)
#'
#' @export
Val.ErrorRates <- function(Network1, Network2, NetworkType = "Association") {
  ## matrices
  matrices_ls <- lapply(list(Network1, Network2), FUN = function(mat) {
    if (NetworkType == "Association") {
      mat[lower.tri(mat)] <- NA
    }
    sign(mat)
  })
  True_mat <- matrices_ls[[1]]
  Inference_mat <- matrices_ls[[2]]

  ## error rates
  errorrates_df <- data.frame(
    Values = c(
      # true positive associations; how many of the inferred positive associations are correctly identified as such? - NaN means no positives in the inferred matrix
      TP = sum((True_mat + Inference_mat) == 2, na.rm = TRUE) /
        sum(Inference_mat == 1, na.rm = TRUE),
      # true negative associations; how many of the inferred negative associations are correctly identified as such? - NaN means no negatives in the inferred matrix
      TN = sum((True_mat + Inference_mat) == -2, na.rm = TRUE) /
        sum(Inference_mat == -1, na.rm = TRUE),
      # true absent; how many of the inferred absent interactions are correctly identified as such? - NaN means no absents in the inferred matrix
      TA = sum((True_mat == 0) + (Inference_mat == 0) == 2, na.rm = TRUE) /
        sum(Inference_mat == 0, na.rm = TRUE),

      # falsely inferred positive; how many of the inferred positive associations are incorrectly identified as such?
      FP = sum((True_mat != 1) + (Inference_mat == 1) == 2, na.rm = TRUE) /
        sum(Inference_mat == 1, na.rm = TRUE),
      # falsely inferred negative; how many of the inferred negative associations are incorrectly identified as such?
      FN = sum((True_mat != -1) + (Inference_mat == -1) == 2, na.rm = TRUE) /
        sum(Inference_mat == -1, na.rm = TRUE),
      # false absent; how many of the inferred absent interactions are correctly identified as such?
      FA = sum((True_mat != 0) + (Inference_mat == 0) == 2, na.rm = TRUE) /
        sum(Inference_mat == 0, na.rm = TRUE),

      # true positive links which were not inferred; how many of the true positive associations were not inferred?
      MP = (sum((True_mat == 1) + (Inference_mat != 1) == 2, na.rm = TRUE) /
        sum(True_mat == 1, na.rm = TRUE)),
      # true negative links which were not inferred; how many of the true negative associations were not inferred?
      MN = (sum((True_mat == -1) + (Inference_mat != -1) == 2, na.rm = TRUE) /
        sum(True_mat == -1, na.rm = TRUE)),
      # true absent links which were not inferred; how many of the true absent associations were not inferred?
      MA = (sum((True_mat == 0) + (Inference_mat != 0) == 2, na.rm = TRUE) /
        sum(True_mat == 0, na.rm = TRUE))
    ),
    Metric = c("TP", "TN", "TA", "FP", "FN", "FA", "MP", "MN", "MA")
  )
  errorrates_df$Values <- errorrates_df$Values * 100
  return(errorrates_df)
}
