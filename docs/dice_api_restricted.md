# Dice Restricted API Reference

This document provides reference documentation for a subset of the public API of the `Dice` library. It is designed to assist agents generating simulation scripts that import the `Dice` solver.

---

## Table of Contents
1. [Core Type Aliases](#core-type-aliases)
2. [Model and Constructors](#model-and-constructors)
3. [File Operations](#file-operations)
4. [State Initialization](#state-initialization)
5. [Solver Execution & Evaluation](#solver-execution--evaluation)
6. [Dynamical Coupling Kernels](#dynamical-coupling-kernels)
7. [Common Usage Examples Explained](#common-usage-examples-explained)

---

## Core Type Aliases

### `SpinConf`
```julia
const SpinConf = Vector{Int8}
```
Represents a discrete spin configuration. Values must be strictly `-1` or `1`.

### `Hybrid`
```julia
const Hybrid = Tuple{SpinConf, Vector{Float64}}
```
Represents the combined state of the network in the relaxed (hybrid) representation. 
- The first element is the discrete spin configuration `SpinConf`.
- The second element is the continuous displacement vector $X \in [-1, 1)^N$.

---

## Model and Constructors

### `Model`
```julia
mutable struct Model
    graph    ::Union{SimpleGraph, SimpleWeightedGraph}
    coupling ::Function
    scale    ::Float64
    silence  ::Int
end
```
Holds the description of the relaxed spin solver.
- `graph`: The network topology.
- `coupling`: The mathematical coupling function.
- `scale`: Magnitude of the simulation time step.
- `silence`: Logging verbosity level (default is `3`, where values $> 3$ show progress/debug info).

#### Constructors
```julia
Model(graph::ModelGraph, coupling::Function)
```
Constructs a model with a default time step scale equal to $1/\Delta(G)$, where $\Delta(G)$ is the maximum vertex degree of the graph.

```julia
Model(graph::ModelGraph, coupling::Function, scale::Float64)
```
Constructs a model with an explicit time step scale.

---

## File Operations

### `load_mtx_graph`
```julia
load_mtx_graph(filename::AbstractString, min_weight::Float64 = 1e-5)::Union{SimpleGraph, SimpleWeightedGraph}
```
Reads a graph from a MatrixMarket (`.mtx`) file. Edges with absolute weights below `min_weight` are ignored. Returns a `SimpleGraph` if all edge weights are identical (within `min_weight` tolerance), or a `SimpleWeightedGraph` otherwise.

### `load_mtx_adjacency`
```julia
load_mtx_adjacency(filename::AbstractString, min_weight::Float64 = 1e-5)::Matrix
```
Convenience function to load a MatrixMarket file and return its full adjacency matrix.

### `save_mtx_graph`
```julia
save_mtx_graph(graph::Union{SimpleGraph, SimpleWeightedGraph}, filename::AbstractString)
```
Saves a graph object to a file in MatrixMarket format.

---

## State Initialization

### `get_random_configuration`
```julia
get_random_configuration(len::Int, p=0.5)::SpinConf
```
Generates a random spin vector of length `len`. Spins are chosen via a Bernoulli process with probability `p` of being `1`, and `1 - p` of being `-1`.

### `get_random_hybrid`
```julia
get_random_hybrid(Nvert::Int, vsize::Float64, p::Float64 = 0.5)::Hybrid
```
Generates a random initial hybrid state `(S, X)` of length `Nvert`.
- The spin component `S` is generated with `get_random_configuration`.
- The continuous displacement component `X` is uniformly distributed in the interval `[-vsize/2, vsize/2]`.

---

## Solver Execution & Evaluation

### `cut`
```julia
cut(graph::Union{SimpleGraph, SimpleWeightedGraph}, conf::SpinConf)::Real
cut(graph::Union{SimpleGraph, SimpleWeightedGraph}, state::Hybrid)::Real
cut(model::Model, state::Hybrid)::Real
```
Evaluates the cut value (the sum of the weights of cut edges) for a given configuration or state. For a hybrid state, the cut is evaluated using only the discrete spin component `S = state[1]`.

### `propagate`
```julia
propagate(model::Model, tmax::Int, start::Hybrid)::Hybrid
```
Advances the simulation of the solver `model` starting from `start` for `tmax` time steps. Returns only the final hybrid state `(S, X)`.

### `trajectories`
```julia
trajectories(model::Model, tmax::Int, start::Hybrid)::Vector{Hybrid}
```
Advances the simulation of the solver `model` starting from `start` for `tmax` time steps. Returns the full history of intermediate hybrid states as a vector of `Hybrid` configurations.

### `propagate_pinned`
```julia
propagate_pinned(model::Model, tmax::Int, start::Hybrid, pinned::Vector{Int64})::Hybrid
```
Advances the simulation for `tmax` steps while keeping the spin variables at indices listed in `pinned` fixed to their initial values.

### `trajectories_pinned`
```julia
trajectories_pinned(model::Model, tmax::Int, start::Hybrid, pinned::Vector{Int64})::Vector{Hybrid}
```
Similar to `propagate_pinned`, but records and returns all intermediate hybrid states.

---

## Dynamical Coupling Kernels

### `model_2_hybrid_coupling`
```julia
model_2_hybrid_coupling(x1::Float64)::Float64
model_2_hybrid_coupling(x1::Float64, x2::Float64)::Float64
```
The sign-based coupling kernel for Model II solvers.
- Single argument: returns `sign(x1)`.
- Two arguments: returns `sign(x1 - x2)`.

---

## Common Usage Examples Explained

### Example 1: Creating a solver for a custom circuit
```julia
function get_onebit_model(parameters) ::Dice.Model
    return Dice.Model(RCAdder.generate_rcadder_graph(1),
                      Dice.model_2_hybrid_coupling,
                      parameters["delta_t"])
end
```
**Explanation**:
1. `RCAdder.generate_rcadder_graph(1)` constructs the graph representing a 1-bit ripple-carry adder circuit.
2. `Dice.model_2_hybrid_coupling` specifies that the simulation will use the sign-based (Model II) dynamics.
3. `parameters["delta_t"]` provides the explicit integration time step scale.
4. Returns a `Dice.Model` ready for propagation.

### Example 2: Agitation and Pinned Propagation Loops
```julia
for _ in 1:parameters["num_agitations"]
    state[2][free] = get_agitation(length(free))
    state = Dice.propagate_pinned(model,
                                  parameters["num_steps"],
                                  state,
                                  pinned)
end
```
**Explanation**:
1. `state` is a `Hybrid` state `(S, X)`. `state[2]` accesses the continuous displacement vector $X$.
2. `free` contains the indices of unpinned (dynamic) nodes.
3. `state[2][free] = get_agitation(length(free))` generates a random perturbation vector and adds noise (agitation) to the continuous coordinates of unpinned nodes to help them escape local energy minima.
4. `Dice.propagate_pinned(model, parameters["num_steps"], state, pinned)` runs the solver dynamics for `num_steps` steps while ensuring the variables at the `pinned` indices remain unchanged.
