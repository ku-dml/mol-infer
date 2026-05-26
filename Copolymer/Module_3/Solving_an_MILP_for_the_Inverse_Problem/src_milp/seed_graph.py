# from .seed_graph_2LMM import *
from .seed_graph_2LMM_polymer import *
# from .seed_graph_2LMM_GNN import *
# from .seed_graph_2LCC import *

def create_seed_graph(sg_type, sg):
	# if sg_type == "2LMM":
	# 	new_sg = SeedGraph_2LMM(name=sg)
	if sg_type == "2LMM_polymer":
		new_sg = SeedGraph_2LMM_Polymer(name=sg)
	# elif sg_type == "2LMM_GNN":
	# 	new_sg = SeedGraph_2LMM_GNN(name=sg)
	# elif sg_type == "2LCC":
	# 	new_sg = SeedGraph_2LCC(name=sg)
	else:
		raise ValueError(f"Unknown seed graph type: {sg_type} !")

	return new_sg