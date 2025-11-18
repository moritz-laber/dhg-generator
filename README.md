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
To generate DHGs you need to specify the number of nodes $n$, the average degree $\bar{k}$, the exponent $\gamma$ of the degree distribution, and the inverse temperature $\beta$. For example,

```python
from dhg import dhg_generator
A, (k,x,c) = dhg_generator(n=10_000, kbar=10, gamma=2.75, beta=2.5)
```

generates a DHG with $n=10^4$ nodes, average degree $\bar{k}=10$, exponent $\gamma=2.75$, and inverse temperature $\beta=2.5$. The graph is returned as an adjacency matrix $A$, and the tuple provides information about the expected degree $k$, angular coordinate $x$, and log-chemical potential $c$.

### Advanced
The `dhg_generator` gives you further control over the graph generation with the following keyword arguments:
* `seed` : The seed for the pseudo-random number generators: These are used when `coordinate_placement` or `edge_placement` are set to `random`, and to seed the Linear Congruential Generator used for pseudo-random coordinate placement when `coordinate_placement` is set to `deterministic`.
* `coordinate_placement` : Can be set to `deterministic` or `random` and controls how the hidden coordinates of nodes are determined. The default is `deterministic`.
* `edge_placement` : Can be set to `deterministic` or `random` and controls how edges are placed. The default is `deterministic`.
* `finite_size_cutoff` : By default `np.inf` is used, and no finite size cut-off is applied. If a finite number is provided, expected degrees are drawn from a truncated Pareto distribution.
* `chemical_potential` : Can be set to `thermodynamic`, in which case the value of the chemical potential in the thermodynamic limit is used or to `numeric`, in which case the value of the chemical potential is determined numerically to match the input average degree exactly. The default is `numeric`.