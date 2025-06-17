# Phase 2

Phase 2 is the inverse QSAR/QSPR phase, aiming at solving an MILP formulation to infer a chemical graph described in Section 3.2.

- MILP_for_2L-GNN
	
	  When given a prediction function obtained from [Phase 1](/2LGNN/Phase_1), a set of rules called *specification* $\sigma$ that describes an abstract structure for a chemical graph, and two real numbers $y_l$ and $y_u$, infer **ONE** chemical graph $C$ such that the predicted property value of $C$ is inside the range $\[y_l, y_u  \]$.
