## Solving an MILP for the Inverse Problem

This folder contains the code used to generate a single chemical graph by solving an MILP formulation for the inverse problem associated with the machine-learning model described in Section 2.3.

The code requires the following input files and parameters:

- files describing the weights/parameters of the prediction functions for the two subsets obtained from [Module 2](/Copolymer/Module_2);
- files describing the topological specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the target range $[y_l, y_u]$.

The installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) is also required. The path to CPLEX should be specified in the configuration file before running the code.

For more details on the MILP formulation, topological specification, and related concepts, please refer to [this preprint](https://arxiv.org/abs/2107.02381).

### Usage

```bash
python infer.py CONFIG.txt OUTPUT
```

Here:

- `CONFIG.txt`: the configuration file specifying the detailed settings;
- `OUTPUT`: the prefix for the output files.

See the sample configuration file [`./sample/sample.yaml`](./sample/sample.yaml) for details.

When the code finishes successfully, it will generate the following files:

- `OUTPUT.lp`: the LP file describing the MILP formulation;
- `OUTPUT.sdf`: the SDF file of the inferred chemical graph, generated when the MILP is feasible;
- `OUTPUT_partition.txt`: the partition file of the inferred chemical graph, generated when the MILP is feasible. This file is used in [Module 4](/HPS/Module_4).

If the solver determines that the MILP is infeasible, then no chemical graph satisfying the given constraints/specifications exists.

### Example

```bash
python infer.py ./sample/sample.yaml test
```