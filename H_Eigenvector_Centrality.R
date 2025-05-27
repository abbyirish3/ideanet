library(igraph)
library(mclust)
library(HyperG)
library(Matrix)
library(RSpectra)

######### CALCULATING H-EIGENVECOTR CENTRALITY ###########

banerjee_coeff <- function(size, max_size) {
  
  sum <- sum(sapply(0:size, function(j) {
    ((-1)^j) * choose(size, j) * ((size - j)^max_size)
  }))
  return(sum)
  
}



# General function subset expansion 
get_gen_coef_subset_expansion <- function(edge_values, node_value, r){
  # initialize the variables
  k <- length(edge_values)
  subset_vector <- c(0)
  subset_lengths <- c(0)
  
  # next, generate all possible subsets of edge values and calculat their sums 
  for (i in 1:k){
    current_length <- length(subset_vector)
    for (t in 1:current_length){
      # add edge_values[i] to the subset t and add its length value to the subset_lengths vector
      subset_vector <- c(subset_vector, subset_vector[t] + edge_values[i])
      subset_lengths <- c(subset_lengths, subset_lengths[t] + 1)
      
    }
  }
  
  # Convert subset lengths to alternating signs based on complement size
  subset_lengths <- (-1)^(k - subset_lengths)
  
  
  # Calculate the weighted sum of (node_value + subset_sum)^r
  total <- sum((node_value + subset_vector)^r * subset_lengths)
  
  # Divide by r factorial
  return(total / factorial(r))
  
  
}


# gen_coeff_fft_array 
#Computes the generating funciton coefficient of order r using the Fast Fourier Transform.
gen_coeff_fft_fast_array <- function(edge_without_node, a, node, l, r){
  
  coefs <- c(1)
  
  # Calculate coefficients for the node being processed
  for (i in 1:(r-1)) {
    coefs <- c(coefs, coefs[length(coefs)] * a[node] / i)
  }
  # Process each node in the edge (except the current node)
  for (u in edge_without_node) {
    temp_coefs <- c(1)
    for (i in 1:(r-l+1)) {
      temp_coefs <- c(temp_coefs, temp_coefs[length(temp_coefs)] * a[u] / i)
    }
    temp_coefs[1] <- 0  # Set first coefficient to 0
    
    # Convolve the coefficient arrays
    coefs <- convolve_r(coefs, temp_coefs)[1:r]
  }
  
  # Return the generating function coefficient
  return(coefs[r])
  
}


convolve_scipy_way <- function(coefs, temp_coefs, r){
  
  coefs <- as.array(coefs)
  temp_coefs <- as.array(temp_coefs)
  
  temp_coefs_rev <- rev(temp_coefs)
  return(convolve(coefs, temp_coefs_rev, type = "open"))
  
  
}


#### Computes the tensor times same vector in all modes but 1.###### 

ttsv1 <- function(node_dict, edge_dict, r, a){
  
  n <- length(node_dict)
  s <- rep(0, n)
  r_minus_1_factorial <- factorial(r-1)
  
  
  for (node in names(node_dict)) {
    node <- as.numeric(node)  # Convert from string to numeric if needed
    c <- 0
    # get all of the edges for the given node
    edges <- node_dict[[as.character(node)]]
    
    for (e in edges){
      
      l <- length(edge_dict[[as.character(e)]])
      alpha <- banerjee_coeff(l, r)
      print("banerjee coeff")
      print(alpha)
      nodes_in_edge <- edge_dict[[as.character(e)]]
      edge_without_node <- nodes_in_edge[nodes_in_edge != node]
      
      if (l == r){
        gen_fun_coef <-prod(a[edge_without_node])
      }
      else if (2^(l - 1) < r * (l - 1)){
        gen_fun_coef <- get_gen_coef_subset_expansion(a[edge_without_node], a[node], r-1)
      }
      else{
        gen_fun_coef <- gen_coeff_fft_fast_array(edge_without_node, a, node, l, r)
        print(gen_fun_coef)
      }
      c <- c + r_minus_1_factorial * l * gen_fun_coef / alpha
      print(c)
    }
    s[node] <- c
    print(s[node])
  }
  
  return(s)
  
}


h_eigenvector_centrality <- function(g, tol=1e-6, max_iter=100){
  
  # input: graph g, tolerance defaulted to 1e-6, max_iter defaulted to 100 
  # output: centrality scores for each node, 1-normalized 
  
  # Extract the two node sets from the bipartite graph
  edges <- V(g)[type==TRUE]
  nodes <- V(g)[type==FALSE]
  
  # Get the actual node IDs/names
  node_ids <- names(nodes)
  edge_ids <- names(edges)
  
  # Check if there aren't any nodes, return an empty list
  if (length(nodes) == 0) {
    return(list())
  }
  
  # Build node_dict: maps each node to the hyperedges it belongs to
  node_dict <- list()
  for (node in nodes) {
    # Get all neighbors of this node (which are hyperedges)
    neighbors <- neighbors(g, node)
    node_dict[[as.character(node)]] <- neighbors
  }
  
  # Build edge_dict: maps each hyperedge to the nodes it contains
  edge_dict <- list()
  for (edge in edges) {
    # Get all neighbors of this hyperedge (which are nodes)
    neighbors <- neighbors(g, edge)
    edge_dict[[as.character(edge)]] <- neighbors
  }
  
  # confirm the graph is not empty 
  if (length(node_ids) == 0) {
    return(list())
  }
  
  # confirm that the graph is connected 
  if (!is_connected(g)) {
    result <- list()
    for (n in node_ids) {
      result[[as.character(n)]] <- NaN
    }
    return(result)
  }
  
  # Determine maximum hyperedge size
  r <- max(sapply(edge_dict, length))
  # Initialize random vector
  num_nodes <- length(node_ids)
  x <- runif(num_nodes)
  # then L-1 normalize
  #x <- x / sum(abs(x)) 
  # normalize
  norm_factor <- sum(abs(x))
  if (norm_factor > 1e-10) {
    x <- x / norm_factor
  } else {
    x <- rep(1/length(x), length(x))  # Fallback to uniform if sum is too small
  }
  
  
  # Create a mapping from node IDs to indices
  node_to_index <- setNames(1:length(nodes), names(nodes))
  
  y <- abs(ttsv1(node_dict, edge_dict, r, x))
  y_scaled <- y^(1/(r-1))
  
  
  converged <- FALSE
  it <- 0
  while (it < max_iter && converged == FALSE){
    
    # Update x
    y_scaled <- y^(1/(r-1))
    x <- y_scaled / sum(abs(y_scaled))  # L1 normalization
    y <- abs(ttsv1(node_dict, edge_dict, r, x))
    
    # now check convergence
    s <- y / (x^(r-1))
    if ((max(s) - min(s)) / min(s) < tol) {
      break
    }
    
    it <- it + 1
    
  }
  
  if (it >= max_iter) {
    warning("Iteration did not converge!")
  }
  
  # Normalize final result
  x <- x / sum(abs(x))
  
  # Create result dictionary with node labels
  result <- setNames(as.list(x), names(nodes))
  
  return(result)
  
}


########## helper functions above ##########
# create graph
g2 <- import_XGI_to_iGraph()
# Plot the bipartite graph
#plot(g2, layout = layout.bipartite, vertex.label = V(g2)$name, 
# main = "Bipartite Graph", vertex.size = 2, edge.width = .1)

# CALCULATE H EIGEN CENTRALITY 


h_centrality2 <- h_eigenvector_centrality(g2)
print(h_centrality2)


convolve_r <- function(in1, in2) {
  # Convert inputs to vectors if they aren't already
  in1 <- as.numeric(in1)
  in2 <- as.numeric(in2)
  
  # Check for scalar inputs
  if (length(in1) == 1 && length(in2) == 1) {
    return(in1 * in2)
  }
  
  # Direct convolution implementation for better precision
  direct_conv <- function(a, b) {
    n_a <- length(a)
    n_b <- length(b)
    n_out <- n_a + n_b - 1
    
    result <- numeric(n_out)
    
    for (i in 1:n_out) {
      sum_val <- 0
      for (j in max(1, i - n_b + 1):min(i, n_a)) {
        k <- i - j + 1
        if (k >= 1 && k <= n_b) {
          sum_val <- sum_val + a[j] * b[k]
        }
      }
      result[i] <- sum_val
    }
    
    return(result)
  }
  
  # For these specific arrays with very small values, direct convolution is likely more accurate
  return(direct_conv(in1, in2))
}


# TESTING 


# first import the graph 
IMPORT_BIPARTITE_EDGELIST_FILEPATH <- "~/Desktop/Mucha/bipartite_edges.csv"

EXPORT_INCIDENCE_MATRIX_FILEPATH <- "~/Desktop/Mucha/bipartite_incidence_matrix.csv"

EXPORT_HEC_DF_FILEPATH <- "~/Desktop/Mucha/hec_scores_r.csv"

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

####### export CEC scores from R #############
hec_df <- data.frame(
  node = names(h_centrality2),
  hec_r = h_centrality_scores
)




