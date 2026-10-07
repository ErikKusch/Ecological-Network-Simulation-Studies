#' Network Inference Accuracy Validation
#'
#' Calculate accuracy of network inference by tallying between two networks the number of similar links by their sign (i.e., absent, positive, and negative links) and divide the number of shared links by the number of all possible links. Networks must be of the same size and contain the same node sets.
#'
#' @param Network1 matrix object. First network in comparison.
#' @param Network2 matrix object. Second network in comparison.
#' @param NetworkType Character. Whether to create an interaction or association matrix. In an interaction matrix, columns act on rows. In an association matrix, entries are mirrored around the diagonal. Possible values: "Association" or "Interaction".
#'
#' @return Percentage of shared links between the two networks out of all possible links.
#'
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#'
#' @examples
#' data("Network_mat")
#' data("Inferred_mat")
#' Val.Accuracy(Network1 = Network_mat, Network2 = Inferred_mat)
#'
#' @export
Val.Accuracy <- function(Network1, Network2, NetworkType = "Association") {
  matrices_ls <- lapply(list(Network1, Network2), FUN = function(mat) {
    if (NetworkType == "Association") {
      mat[lower.tri(mat)] <- NA
    }
    sign(mat)
  })
  eq <- matrices_ls[[1]] == matrices_ls[[2]]
  sum(eq, na.rm = TRUE) / sum(!is.na(eq)) * 100
}
