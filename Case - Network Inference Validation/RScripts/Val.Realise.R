#' Reduction of Network to Realisable Component
#'
#' Takes an association or interaction network object as well as the bioclimatic niche preferences of each species to compute which links stored in the network can be realised. The cutoff is calculated as the difference between the niche preferences of species pairs vs. Effect_Dis + Env_sd.
#'
#' @param Network_mat An matrix object with association/interaction strength stored as cell values. Output of Sim.Network().
#' @param NetworkType Character. Whether to create an interaction or association matrix. In an interaction matrix, columns act on rows. In an association matrix, entries are mirrored around the diagonal. Possible values: "Association" or "Interaction".
#' @param Trait_means Numeric. Mean for bioclimatic niche preference of each species.
#' @param Effect_Dis Numeric. Distance within which neighbouring individuals interact with a focal individual.
#' @param Env_sd Numeric. Habitat suitability in death rate function. Higher values allow individuals to persist in areas of greater environmental maladaptation.
#'
#' @return A matrix object with association/interaction strength stored as cell values.
#'
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#'
#'
#' @examples
#' data("Niches_vec")
#' data("Network_mat")
#' Val.Realise(Network_mat = Network_mat, Trait_means = Niches_vec)
#'
#' @export
Val.Realise <- function(Network_mat,
                        NetworkType = "Association",
                        Trait_means,
                        Effect_Dis = 0.5,
                        Env_sd = 2.5) {
  ## Vertex names
  # if (is.null(igraph::V(Network_igraph)$names)) {
  #   igraph::V(Network_igraph)$names <- paste0("Sp_", igraph::V(Network_igraph))
  # }

  ## Mode of network
  # mode <- ifelse(is_directed(Network_igraph), "directed", "undirected")

  ## Network Adjacency matrix
  # Network_Realised <- as.matrix(as_adjacency_matrix(Network_igraph, attr = "weight"))
  # colnames(Network_Realised) <- rownames(Network_Realised) <- igraph::V(Network_igraph)$names

  ## Trait Adjacency Matrix
  # Figuring out trait differences between potentially interacting species
  SPTrait_df <- data.frame(Trait_means)
  SPTrait_df$Species <- rownames(SPTrait_df)
  colnames(SPTrait_df) <- c("Trait", "Species")

  SPTrait_df <- SPTrait_df[match(colnames(Network_mat), SPTrait_df$Species), ]
  SPTrait_mat <- abs(outer(SPTrait_df$Trait, SPTrait_df$Trait, "-"))
  colnames(SPTrait_mat) <- rownames(SPTrait_mat) <- colnames(Network_mat)

  TraitDiff_mat <- SPTrait_mat

  ## limitting to realised interactions
  Network_mat[(TraitDiff_mat) > (Env_sd + Effect_Dis)] <- 0 # anything greater apart in enviro pref than the interaction window + environmental sd cannot be mat

  ## make into matrix
  Real_mat <- Network_mat
  diag(Real_mat) <- NA
  if (NetworkType == "Association") {
    Real_mat[lower.tri(Real_mat)] <- NA
  }

  ## make into igraph object
  # Network_Real <- graph_from_adjacency_matrix(
  #   adjmatrix = Network_Realised,
  #   mode = mode,
  #   weighted = TRUE,
  #   diag = FALSE
  # )

  return(Real_mat)
}
