# Module 3

Module 3 consists of the stage of solving an MILP formulation for the inverse problem to infer a chemical graph described in Section 6.2, paragraphs "Solving an MILP for the Inverse Problem" and "Generating Neighbor Solutions".

- Solving an MILP for the Inverse Problem
  
  When given a prediction function obtained from [Module 2](/HPS/Module_2), a set of rules called *topological configuration* $\sigma$ that describes an abstract structure for a chemical graph, and two real numbers $y_l$ and $y_u$, infer **ONE** chemical graph $C$ such that the predicted property value of $C$ is inside the range $\[y_l, y_u  \]$.

- Generating Neighbor Solutions

  We add a set of new linear constraints to the MILP above in order to construct more solutions systematically.
