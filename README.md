# DHG Generator
This repository contains code for the **Deterministic Hyperbolic Graph (DHG) Generator** from the paper *Deterministic Construction of Typical Networks in Network Models* and a brief tutorial on how to use it.

## Installation 
To use the DHG generator you can directly `pip install` from this repository using

```bash
pip install "git+https://github.com/moritz-laber/dhg-generator"
```

If you want to experiment with the tutorial, you might want to clone the repository and install additional dependencies via

```bash
git clone https://github.com/moritz-laber/dhg-generator
cd dhg-generator
pip install -e ".[tutorial]"
```

## Usage

### Basic
The Deterministic Hyperbolic Graph (DHG) generator `dhg_generator` has four required inputs

* `n` - the number of nodes $n$,
* `kbar` - the average degree $\bar{k}$,
* `gamma` - the exponent $\gamma$ of the power law tail of the degree distribution,
* `beta` - the inverse temperature $\beta$,

and returns

* `A` - the adjacency matrix $A$ of the generated graph in `scipy.sparse.coo_matrix` format.
* a tuple `(kappa, x, mu)` of further properties:
    * `kappa` - an `numpy.ndarray` such that entry $i$ is the expected degree $\kappa_i$ of node $i$,
    * `x` - an `numpy.ndarray` such that entry $i$ is the angular coordinate $x_i\in[0,n)$ of node $i$,
    * `mu` - the chemical potential $\mu$ of the graph.

To generate a graph you can use the following function call:

```python
from dhg import dhg_generator
A, (kappa, x, mu) = dhg_generator(n=10_000, kbar=10, gamma=2.75, beta=2.5)
```

This generates a DHG with $n=10^4$ nodes, average degree $\bar{k}=10$, exponent $\gamma=2.75$, and inverse temperature $\beta=2.5$. The details of the graph construction are explained in Appendix G of the paper *Deterministic Construction of Typical Networks in Network Models*. The outputs `kappa`, `x`, and `mu` can be used to compute further properties of the graph that are not accessible from the adjacency matrix alone, e.g., the energy of the graph.

### Advanced
The `dhg_generator` gives you further control over the graph generation with the following keyword arguments:
* `seed` : The seed for the pseudo-random number generators: These are used when `coordinate_placement` or `edge_placement` are set to `random`, and to seed the Linear Congruential Generator used for pseudo-random coordinate placement when `coordinate_placement` is set to `deterministic`.
* `coordinate_placement` : Can be set to `deterministic` or `random` and controls how the hidden coordinates of nodes are determined. The default is `deterministic`.
* `edge_placement` : Can be set to `deterministic` or `random` and controls how edges are placed. The default is `deterministic`.
* `finite_size_cutoff` : By default `np.inf` is used, and no finite size cut-off is applied. If a finite number is provided, expected degrees are drawn from a truncated Pareto distribution.
* `chemical_potential` : Can be set to `thermodynamic`, in which case the value of the chemical potential in the thermodynamic limit is used or to `numeric`, in which case the value of the chemical potential is determined numerically to match the input average degree exactly. The default is `numeric`.

### Using `dhg_generator` with network analysis packages

From the adjacency matrix, you can easily construct graph objects in popular network analysis packages such as [graph-tool](https://graph-tool.skewed.de/)

```python
from dhg import dhg_generator
import graph_tool as gt
from graph_tool import *
A, _ = dhg_generator(n=10_000, kbar=10, gamma=2.75, beta=2.5)
g = gt.Graph(directed=False)
g.add_nodes(A.shape[0])
g.add_edge_list(zip(A.row, A.col))
gt.generation.remove_parallel_edges(g)
```

or [NetworkX](https://networkx.org/documentation/stable/index.html)

```python
from dhg import dhg_generator
A, _ = dhg_generator(n=10_000, kbar=10, gamma=2.75, beta=2.5)
g = nx.from_scipy_sparse_array(A)
```