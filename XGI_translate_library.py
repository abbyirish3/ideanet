'''
Todd Rosenbaum, 26 May 2025
Ideanet research
Functions for import/exprt with XGI/iGraph
'''

import xgi
import matplotlib.pyplot as plt
import pandas as pd

##### GLOBAL VARIABLES #####
IMPORT_INCIDENCE_MATRIX_FILEPATH = "/Users/toddrosenbaum/Desktop/Mucha/bipartite_incidence_matrix.csv"
EXPORT_BIPARTITE_EDGELIST_FILEPATH = "/Users/toddrosenbaum/Desktop/Mucha/bipartite_edges.csv"

'''
Import from a bipartite incidence matrix and represent as XGI hypergraph
Input: filepath to bipartite incidence matrix
Return: XGI hypergraph
'''
def import_iGraph_to_XGI():
    # Create incidence matrix
    incidence_matrix = pd.read_csv(IMPORT_INCIDENCE_MATRIX_FILEPATH, index_col=0)
    # Convert the incidence matrix to a dictionary format suitable for xgi
    # Rows represent one node type, columns represent the other
    hypergraph_dict = {
        col: incidence_matrix.index[incidence_matrix[col] == 1].tolist()
        for col in incidence_matrix.columns
    }
    # Create and return hypergraph
    return xgi.Hypergraph(hypergraph_dict)


'''
Export from an XGI hypergraph to a bipartite edgelist
Input: XGI hypergraph
Return: none
Do: create edgelist file
'''
def export_XGI_to_iGraph(hypergraph):
    # Create edgelist
    edge_list = xgi.to_bipartite_edgelist(hypergraph)
    # save bipartite edge_list as csv
    with open(EXPORT_BIPARTITE_EDGELIST_FILEPATH, "w") as f:
        f.write("node,hyperedge\n")

        for edge in edge_list:
            f.write(f"{edge[0]},{edge[1]}\n")


###########################################
########## Helper function above ##########
########## execute code below #############
###########################################

########## Test hypergraph ##########
'''
H = xgi.Hypergraph()
H.add_nodes_from([
    (1, {"label": "1"}),
    (2, {"label": "2"}),
    (3, {"label": "3"}),
    (4, {"label": "4"}),
    (5, {"label": "5"}),
])
H.add_edges_from([
    ({1, 2, 3}, {"weight": 1.0}),
    ({2, 3, 4}, {"weight": 1.5}),
    ({1, 4, 5}, {"weight": 2.0}),
    ({3, 5}, {"weight": 0.5}),
])
'''

########## XGI hypergraph ##########
# I have found the following data sets to be the most helpful: plant-pollinator-mpl-014, senate-bills, senate-committees
# see https://xgi.readthedocs.io/en/stable/xgi-data.html for more
H = xgi.load_xgi_data(dataset="senate-committees")

# export H
export_XGI_to_iGraph(H)
# import J
J = import_iGraph_to_XGI()


########## XGI clustering coefficients ##########
print(xgi.clustering_coefficient(H))
print(xgi.local_clustering_coefficient(H))
print(xgi.two_node_clustering_coefficient(H, "max")) # type = (max, min, union)


########## hypergraph visualization ##########
# plots the hypergraph H and re-imported hypergraph J side by side
# Create subplots for side-by-side display
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# Display original hypergraph H
axs[0].set_title("Original Hypergraph (H)")
xgi.draw(H, ax=axs[0])
# Display re-imported hypergraph J
axs[1].set_title("Re-imported Hypergraph (J)")
xgi.draw(J, ax=axs[1])
# Show both drawings together
plt.show()
