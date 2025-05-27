library(igraph)
library(mclust)
library(HyperG)
library(Matrix)
library(RSpectra)

# first import the graph 
IMPORT_BIPARTITE_EDGELIST_FILEPATH <- "~/Desktop/bipartite_edges.csv"
EXPORT_INCIDENCE_MATRIX_FILEPATH <- "~/Desktop/bipartite_incidence_matrix.csv"

# read from edgelist and create/plot bipartite graph
import_XGI_to_iGraph <- function() {
  # read edgelist
  edges = read.csv(IMPORT_BIPARTITE_EDGELIST_FILEPATH, header = TRUE)
  # modify hyperedges in the table to avoid igraph considering hyperedge same as node
  edges$hyperedge <- paste0("T", edges$hyperedge)
  # Create the bipartite graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  # Set the bipartite vertex type (FALSE for 'from' nodes, TRUE for 'to' nodes)
  V(g)$type <- c(rep(FALSE, length(unique(edges$node))), rep(TRUE, length(unique(edges$hyperedge))))
  
  return(g)
}

# write incidence matrix of bipartite graph
export_iGraph_to_XGI <- function(bipartite_graph) {
  # Convert to incidence matrix
  inc_matrix <- as_biadjacency_matrix(bipartite_graph, types=)
  # Save as CSV
  write.csv(as.matrix(inc_matrix), EXPORT_INCIDENCE_MATRIX_FILEPATH, row.names = TRUE)
}


########## helper functions above ##########

# create graph
g <- import_XGI_to_iGraph()
# Plot the bipartite graph
plot(g, layout = layout.bipartite, vertex.label = V(g)$name, 
     main = "Bipartite Graph", vertex.size = 2, edge.width = .1)
# export_iGraph_to_XGI
export_iGraph_to_XGI(g)


######## compute the clique-eigenvector centrality #########

# get biadjacency matrix: rows are nodes from one sie of bipartite, columns are edges: Bij = 1 if node i is in hyperedge j
types <- V(g)$type
B <- as_biadjacency_matrix(g, types = types, sparse = TRUE) # B will be |U| x |V| matrix

# project onto one side: W = BB^T
W <- B %*% t(B)  # |U| x |U| symmetric matrix, Wuv = # hypereges shared by nodes u and v

# compute largest eigenvector of W
cec <- eigs_sym(W, k = 1, which = "LA")  # LA returns largest real eigenvalue
centrality_scores <- cec$vectors[, 1]
centrality_scores <- centrality_scores / sum(centrality_scores)  # normalize
names(centrality_scores) <- rownames(W)

print(centrality_scores)
print(max(centrality_scores))


####### export CEC scores from R #############
cec_df <- data.frame(
  node = names(centrality_scores),
  cec_r = centrality_scores
)
write.csv(cec_df, "~/Desktop/cec_scores_r.csv", row.names = FALSE)


####### Compute ZEC #############
# build tensor representation from incidence matrix
# extract_tensor_entries <- function(B, m) {
#   tensor_entries <- list()
#   edge_indices <- colnames(B)
#   idx <- 1
#   
#   for (e in 1:ncol(B)) {
#     nodes_in_edge <- which(B[, e] != 0)
#     if (length(nodes_in_edge) < m) {
#       warning(paste("Edge", e, "has fewer than", m, "nodes"))
#       next
#     }
#     # all combinations of m nodes in this hyperedge
#     combos <- combn(nodes_in_edge, m)
#     for (col in 1:ncol(combos)) {
#       indices <- as.integer(combos[, col])
#       entry <- list(indices = indices, val = 1)
#       tensor_entries[[idx]] <- entry
#       idx <- idx + 1
#     }
#   }
#   return(tensor_entries)
# }
extract_tensor_entries <- function(B, m) {
  tensor_entries <- list()
  idx <- 1
  
  for (e in 1:ncol(B)) {
    nodes_in_edge <- which(B[, e] != 0)
    r_e <- length(nodes_in_edge)
    
    if (r_e < m) {
      warning(paste("Edge", e, "has fewer than", m, "nodes"))
      next
    }
    
    combos <- combn(nodes_in_edge, m)
    alpha_e <- choose(r_e - 1, m - 1)  # this is α(e)
    
    for (col in 1:ncol(combos)) {
      indices <- as.integer(combos[, col])
      entry <- list(indices = indices, val = 1 / alpha_e)
      tensor_entries[[idx]] <- entry
      idx <- idx + 1
    }
  }
  
  return(tensor_entries)
}





# tensor multiplication function
# tensor_mult <- function(tensor_entries, c, m) {
#   n <- length(c)
#   result <- rep(0, n)
#   
#   for (entry in tensor_entries) {
#     indices <- entry$indices
#     val <- entry$val
#     
#     i1 <- indices[1]
#     other_indices <- indices[-1]
#     prod_c <- prod(c[other_indices])
#     
#     result[i1] <- result[i1] + val * prod_c
#   }
#   
#   return(result)
# }
tensor_mult <- function(tensor_entries, c, m) {
  n <- length(c)
  result <- rep(0, n)
  
  for (entry in tensor_entries) {
    indices <- entry$indices
    val <- entry$val
    prod_c <- prod(c[indices])
    
    for (i in indices) {
      result[i] <- result[i] + val * prod_c / m  # distribute equally
    }
  }
  return(result)
}

# power iteration
z_eigenvector <- function(tensor_entries, m, max_iter=100, tol=1e-6) {
  n <- max(sapply(tensor_entries, function(e) max(e$indices)))
  c <- rep(1 / n, n)  # uniform init with ℓ1 norm
  
  for (iter in 1:max_iter) {
    new_c <- tensor_mult(tensor_entries, c, m)
    new_c <- new_c / sum(abs(new_c))  # switch to l1 norm
    
    denom <- min(abs(c))
    if (denom < 1e-10) {
      denom <- 1e-10
    }
    if ((max(abs(new_c - c)) / denom) < tol) {
      cat("Converged at iteration", iter, "\n")
      break
    }
    }
    c <- new_c
  return(c / sum(abs(c)))  # l1 normalization
}


## example using matrix B from earlier
# Extract tensor entries
m <- 3
tensor_entries <- extract_tensor_entries(B, m)

# Compute ZEC
zec <- z_eigenvector(tensor_entries, m)

# Label results
names(zec) <- rownames(B)
print(zec)



