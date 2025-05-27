# Todd Rosenbaum, 26 May 2025
# Ideanet research
# Functions for import/exprt with XGI/iGraph

library(igraph)

########## Global variables ##########
IMPORT_BIPARTITE_EDGELIST_FILEPATH <- "/Users/toddrosenbaum/Desktop/Mucha/bipartite_edges.csv"
EXPORT_INCIDENCE_MATRIX_FILEPATH <- "/Users/toddrosenbaum/Desktop/Mucha/bipartite_incidence_matrix.csv"


# read from edgelist and create bipartite graph
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
  inc_matrix <- as_biadjacency_matrix(bipartite_graph)
  # Save as CSV
  write.csv(as.matrix(inc_matrix), EXPORT_INCIDENCE_MATRIX_FILEPATH, row.names = TRUE)
}


# one-mode projection
bipartite_to_unipartite <- function(bipartite_graph) {
  # create projection
  g_proj <- bipartite_projection(bipartite_graph)
  # assign which are nodes and which are edges
  g_unipartite <- g_proj$proj1
  
  return(g_unipartite)
}


# clustering coefficient
clustering_coefficient <- function(unipartite_graph){
  return(transitivity(g_unipartite, type = "local", isolates = "zero"))
}

###############################################################################################################
########## THIS WORKS FOR SMALL TOY-GRAPHS BUT TAKES EXTREMELY LONG FOR MORE COMPLICATED HYPERGRAPHS ##########
###############################################################################################################
# local clustering coefficient
local_clustering_coefficient <- function(g) {
  # Identify original nodes and hyperedges
  orig_idxs <- which(V(g)$type == FALSE)
  edge_idxs <- which(V(g)$type == TRUE)
  orig_names <- V(g)$name[orig_idxs]
  edge_names <- V(g)$name[edge_idxs]
  
  # Precompute name-to-id map
  name_to_vid <- setNames(seq_along(V(g)), V(g)$name)
  
  # Precompute neighbors for all nodes
  neighbor_map <- lapply(seq_along(V(g)), function(v) neighbors(g, v))
  names(neighbor_map) <- V(g)$name
  
  # Precompute edge members (edge name → original nodes)
  edge_members <- lapply(edge_idxs, function(ei) {
    vids <- neighbor_map[[V(g)$name[ei]]]
    vids[V(g)[vids]$type == FALSE]
  })
  names(edge_members) <- edge_names
  
  # Precompute node memberships (node name → incident edge names)
  node_memberships <- lapply(orig_idxs, function(vi) {
    vids <- neighbor_map[[V(g)$name[vi]]]
    vids[V(g)[vids]$type == TRUE]
  })
  names(node_memberships) <- orig_names
  
  result <- setNames(numeric(length(orig_names)), orig_names)
  
  for (v_name in orig_names) {
    evs <- node_memberships[[v_name]]
    dv <- length(evs)
    if (dv <= 1) {
      result[[v_name]] <- 0
      next
    }
    
    ev_names <- V(g)$name[evs]
    pairs <- combn(ev_names, 2, simplify = FALSE)
    total_eo <- 0
    
    for (pr in pairs) {
      e1v <- pr[1]
      e2v <- pr[2]
      
      D1 <- setdiff(edge_members[[e1v]], edge_members[[e2v]])
      D2 <- setdiff(edge_members[[e2v]], edge_members[[e1v]])
      U <- union(D1, D2)
      eo <- 0
      
      if (length(U) > 0) {
        # Get neighbor sets
        neighD1 <- unique(unlist(lapply(D1, function(nm) {
          hvs <- neighbor_map[[nm]]
          ov <- unlist(lapply(hvs, function(hv) {
            nbrs <- neighbor_map[[V(g)$name[hv]]]
            V(g)$name[nbrs[V(g)[nbrs]$type == FALSE]]
          }))
          ov
        })))
        
        neighD2 <- unique(unlist(lapply(D2, function(nm) {
          hvs <- neighbor_map[[nm]]
          ov <- unlist(lapply(hvs, function(hv) {
            nbrs <- neighbor_map[[V(g)$name[hv]]]
            V(g)$name[nbrs[V(g)[nbrs]$type == FALSE]]
          }))
          ov
        })))
        
        eo <- (length(intersect(neighD1, D2)) + length(intersect(neighD2, D1))) / length(U)
      }
      
      total_eo <- total_eo + eo
    }
    
    result[[v_name]] <- 2 * total_eo / (dv * (dv - 1))
  }
  
  return(result)
}


# two_node clustering coefficient
two_node_clustering_coefficient <- function(g, kind = c("union","min","max")) {
  kind      <- match.arg(kind)
  originals <- V(g)[V(g)$type == FALSE]
  result    <- numeric(length(originals))
  names(result) <- originals$name
  
  for (i in seq_along(originals)) {
    u <- originals[i]
    
    # 1) hyperedges incident on u
    H_u <- neighbors(g, u)[ V(g)[neighbors(g,u)]$type ]
    # 2) original‐mode neighbors via those hyperedges
    via_h  <- unlist(lapply(H_u, function(e) neighbors(g, e)))
    neighs <- setdiff(unique(via_h[ V(g)[via_h]$type == FALSE ]), u)
    
    if (length(neighs)==0) {
      result[i] <- 0
    } else {
      sims <- vapply(neighs, function(v) {
        H_v <- neighbors(g, v)[ V(g)[neighbors(g,v)]$type ]
        uv_cc(H_u, H_v, kind)
      }, numeric(1))
      result[i] <- mean(sims)
    }
  }
  
  return(result)
}

# helper for two_node clustering coefficient
uv_cc <- function(H_u, H_v, kind = c("union","min","max")) {
  kind <- match.arg(kind)
  num   <- length(intersect(H_u, H_v))
  denom <- switch(kind,
                  union = length(union(H_u, H_v)),
                  min   = min(length(H_u), length(H_v)),
                  max   = max(length(H_u), length(H_v))
  )
  if (denom == 0) return(NaN)
  num/denom
}