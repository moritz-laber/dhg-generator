import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from typing import Tuple

def lcg(n: int, x0: int = 42, a: int = 7**5 , m: int = 2**31 - 1) -> np.ndarray:
    """Create n pseudo-random numbers using a Lehmer generator.
    
    Input
    n  - number of random numbers.
    x0 - seed.
    a  - multiplier.
    m  - modulus.
    
    Output
    numbers - array of n pseudo random numbers.
    """

    # set initial condition
    x = x0

    # generate pseudo-random numbers
    numbers = []
    for _ in range(n):

        # update the state
        x = (a*x) % m

        numbers.append(x)
    
    return np.asarray(numbers)


def pareto_inverse_cdf(u: np.ndarray, xbar: float, gamma: float, n:int=np.inf) -> np.ndarray:
    """Inverse cumulative distribution function of the Pareto distribution.
    
    Input
    u - values at which to evaluate the inverse cdf.
    xbar - the mean of the Pareto distribution.
    gamma - the exponent of the probability density function.
    n - the number of nodes for the Pareto distribution with finite size cutoff. Set np.inf for no cutoff.

    Output
    x - the values of the inverse cdf. If u is uniformly distributed, x is Pareto distributed.
    """

    if n < np.inf:
        x = xbar * ((gamma - 2) / (gamma - 1))*((1-1/n)/(1-n**((gamma-2)/(1-gamma))))*(1 - (1 - 1/n)*u)**(1 / (1 - gamma))
    else:
        x = xbar * ((gamma - 2) / (gamma - 1)) * (1 - u)**(1 / (1 - gamma))

    return x

def exponential_inverse_cdf(u: np.ndarray, lam: float) -> np.ndarray:
    """Inverse cumulative distribution function of the exponential distribution.
    
    Input
    u - values at which to evaluate the inverse cdf.
    lam - the rate parameter of the exponential distribution.

    Output
    x - the values of the inverse cdf. If u is uniformly distributed, x is exponentially distributed.
    """

    lam = float(lam)
    
    return (-1.0/lam) * np.log(1.0 - u)

def deterministic_coordinates(n:int, z:float=0.5, x0:int=42) -> Tuple[np.ndarray]:
    """Create n pseudo-random points in [0, 1] x [0, 1].
    
    Input
    n - number of points for which to create coordinates
    z - shift of the lattice points
    x0 - the seed for the pseudo-random permutations
    
    Output
    x,y - arrays of coordinates
    """

    # create a lattice that uniformly covers [0, 1]
    u = np.arange(1 - z , n + 1 - z) / n
    d_nn  = exponential_inverse_cdf(u, n)

    # create pseudo-random permutations
    r = lcg(4*n, x0=x0)
    
    sigma1 = np.argsort(r[:n])
    sigma2 = np.argsort(r[n:2*n])
    sigma3 = np.argsort(r[2*n:3*n])
    sigma4 = np.argsort(r[3*n:4*n])

    # create the x and y coordinates from the nearest neighbor distances
    x = np.cumsum(d_nn[sigma1])
    y = np.cumsum(d_nn[sigma2])

    x = x[sigma3]
    y = y[sigma4]

    return x, y

def random_coordinates(n:int, rng:np.random.Generator) -> Tuple[np.ndarray]:
    """Create n random points in [0, 1] x [0, 1].
    
    Input
    n - number of points
    rng - a numpy pseudorandom number generator

    Output
    x,y - coordinates of points.    
    """
    
    # create random coordinates
    x = rng.random(n)
    y = rng.random(n)

    return x, y


def objective_function(c: float, n: int, kbar: float, beta: float, chi_beta: np.ndarray) -> Tuple[float]:
    """Measures the deviation of the average degree from its expected value given a set of 
    coordinates and a log chemical potential.
    
    Input
    c - the log chemical potential
    n - the number of nodes
    kbar - the expected degree
    beta - the inverse temperature
    chi_beta - log scaled hyperbolic distance exponentiated by inverse temperature
    
    Output
    r - residual of the optimization procedure
    dr - derivative of the residual of the optimization procedure   
    """

    # exponentiate for better readability
    c_beta = c**beta

    # calculate the connection probability
    p = 1. / (1. + chi_beta/c_beta)

    # calculate the derivative of the connection probability w.r.t. c
    dpdc = (beta/c) * (chi_beta/c_beta) / (1. + chi_beta/c_beta)**2

    # calculate the residual of the optimization
    r = 2.0 * np.sum(p) / n - kbar

    # calculate the change in residual w.r.t. c 
    drdc = 2.0 * np.sum(dpdc) / n

    return r, drdc

def newton(x0:np.ndarray, objective:callable, atol:float=1e-8, maxiter:int=int(1e5))->Tuple[np.ndarray, bool]:
    """Newton's method for finding roots of multivariate scalar function f(x).
    
    Input
    x0 - the initial guess
    objective - callable that takes x and returns f(x), df/dx
    atol - absolute tolerance for convergence (default: 1e-8)
    maxiter - maximum number of iterations (defualt: 1e5)
    
    Output
    x - the root found
    converged - whether the method converged
    """ 

    # set the initial condition
    x = x0
    converged = False

    for _ in range(int(maxiter)):

        # evaluate function and derivative
        f_x, dfdx_x = objective(x)

        # check convergence
        if np.abs(f_x) < atol:
            converged = True
            break
        
        # take one Newton step
        x = x - f_x / dfdx_x

    if not converged:
        Warning("Newton's method did not converge within the maximum number of iterations. Returning the last iterate.")

    return x, converged

def deterministic_placement(p: np.ndarray) -> np.ndarray:
    """Method for deterministic edge-placement given the probability of each edge.
    
    Input
    p - array of edge probabilities 

    Output
    edges_idx - 1d array of indices of p for which an edge exists 
    """

    # accumulate the probabilities
    p_cummulative = np.cumsum([0, *p])

    # determine the number of links
    mclosest = int(np.round(p_cummulative[-1]))

    # create a grid
    x = np.arange(0, mclosest) + 0.5

    # place the links pseudo-randomly
    edges_idx = np.searchsorted(p_cummulative, x, side='right') - 1

    return edges_idx


def random_placement(p: np.ndarray, rng:np.random.Generator) -> np.ndarray:
    """Method for random placement of edges given probability of each edge.
    
    Input
    p - the edge probabilities

    Output
    edges_idx - 1d array of indices of p for which an edge exists
    """

    # determine the shape of the adjacency matrix
    mmax = p.shape[0]

    # draw n uniform random numbers on [0,1] 
    r = rng.random(mmax)

    # decide whether edges should be present or absent
    edges_idx = np.where(r < p)[0]

    return edges_idx


def dhg_generator(
        n: int,
        kbar: float,
        gamma: float,
        beta: float,
        seed:int=42,
        coordinate_placement:str = "deterministic",
        edge_placement:str = "deterministic",
        finite_size_cutoff:bool = False,
        chemical_potential:str = "thermodynamic",
        ) -> Tuple[coo_matrix, Tuple[np.ndarray, np.ndarray, float]]:
    """Generates a deterministic hyperbolic graph (DHG) with n nodes, average degree kbar, exponent gamma of the degree distribution, and inverse temperature beta.
    
    Input
    n     - the number of nodes 
    kbar  - the average degree
    gamma - the exponent of the degree distribution
    beta  - the inverse temperature
    seed - seed for the pseudorandom number generator
    coordinate_placement - whether to place coordinates deterministically ("deterministic") or randomly ("random")
    edge_placement - whether to place the edges deterministically ("deterministic") or randomly ("random")
    finite_size_cutoff - whether to use the truncated Pareto distribution with finite size cutoff (True) or not (False)
    chemical_potential - whether to use the chemical potential in the thermodynamic limit ("thermodynamic") or the numerically determined one ("numeric")

    Output
    A - adjacency matrix of the graph
    (kappa, x, mu) - additional properties: the expected degree kappa_i (popularity) and angular coordinate x_i (similarity) of each node i as welll as the chemical potential mu. 
    """

    # generate subseeds for the random case
    if coordinate_placement=="random" or edge_placement=="random":

        seed_rng = np.random.default_rng(seed)
        coordinate_seed = seed_rng.integers(low=1, high=2**31 - 1)
        edge_seed = seed_rng.integers(low=1, high=2**31 - 1)

    # create a coordinate assignment
    if coordinate_placement == "deterministic":
        
        x, y = deterministic_coordinates(n, x0=int(seed))
    
    elif coordinate_placement == "random":

        x, y = random_coordinates(n, rng=np.random.default_rng(coordinate_seed))
    
    else:

        raise ValueError("Unknown coordinate_placement. Only 'deterministic' or 'random' supported.")

    # transform coordinates to degree-like latent variables
    if finite_size_cutoff==True:
        kappa = pareto_inverse_cdf(y, kbar, gamma, n=n)
    else:
        kappa = pareto_inverse_cdf(y, kbar, gamma)

    
    # calculate the distance respecting periodic boundary conditions
    idx_i, idx_j = np.triu_indices(n, k=1)
    chi_beta = np.power(n*(0.5 - np.abs(0.5 - np.abs(x[idx_i] - x[idx_j]))) / (kappa[idx_i]*kappa[idx_j]), beta)

    # determine the chemical potential in the thermodynamic limit
    c0 = (np.sin(np.pi/beta)/(np.pi/beta)) / (2*kbar)
    
    # determine the chemical potential with the selected method
    if chemical_potential=="thermodynamic":
        
        c = c0

    elif chemical_potential=="numeric":

        c, converged = newton(x0=c0, objective=lambda c: objective_function(c, n, kbar, beta, chi_beta))

        if converged:
            c

        else:
            c = c0
            raise Warning("Could not determine the chemical potential. Falling back to the value in the value in the thermodynamic limit.")
    
    else:
        raise ValueError("Unknown chemical_potential. Only 'thermodynamic' or 'numeric' supported.")

    # calculate the connection probability
    p = 1. / (1. + chi_beta / c**beta)

    # place the edges
    if edge_placement == "deterministic":

        edges_idx = deterministic_placement(p)
    
    elif edge_placement == "random":

        edges_idx = random_placement(p, rng=np.random.default_rng(edge_seed))

    else:

        raise ValueError("Unknown edge_placement. Only 'deterministic' or 'random' supported.")

    # create the adjacency matrix
    A = coo_matrix((np.ones(edges_idx.shape[0]), (idx_i[edges_idx], idx_j[edges_idx])), shape=(n,n), dtype=np.int32)
    A = A + A.T
    A.tocoo()

    # rescale to S1 model for output
    mu = np.log(c)
    x *= n

    return A, (kappa, x, mu)