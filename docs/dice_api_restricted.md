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
7. [Saving Metadata with Simulation Results](#saving-metadata-with-simulation-results)
8. [Common Usage Examples Explained](#common-usage-examples-explained)

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

## Saving Metadata with Simulation Results

The library supports embedding machine-readable JSON metadata records in simulated data files (such as those generated by `Dice.save_trajectories`). To maintain strict backward compatibility, these metadata records are formatted as JSON objects embedded inside comment lines prefixed with `#@`. Legacy parsers ignore these comments, while enhanced metadata-aware parsers extract and interpret them.

### `Dice.save_trajectories`
```julia
Dice.save_trajectories(io::IO, traj::Vector{Hybrid}, time::AbstractVector{<:Real}, header=[])
Dice.save_trajectories(path::AbstractString, traj::Vector{Hybrid}, time::AbstractVector{<:Real}, header=[])
Dice.save_trajectories(traj::Vector{Hybrid}, header=[])
```
Writes a relaxed-spin trajectory dataset to a stream or file path. If `time` is omitted, it defaults to a `0`-based integer step index sequence.

### Header Construction Helpers

The following helper functions under the `Dice` namespace are provided to build items for the `header` array parameter:

* **`Dice.comment(text::AbstractString)`**
  Generates an ordinary human-readable comment line starting with `# `.
* **`Dice.format_record(; name="relaxed-spins", version=1)`**
  Generates the standard relaxed-spin format version declaration.
* **`Dice.created(value::AbstractString)` / `Dice.created(dt::DateTime)`**
  Generates a creation timestamp. `DateTime` objects are formatted as UTC ISO strings.
* **`Dice.title(text::AbstractString)`**
  Generates a title metadata record.
* **`Dice.notes(text::AbstractString)`**
  Generates general notes/description metadata.
* **`Dice.node_indexing(base::Integer = 1)`**
  Declares whether node indexing uses `0`-based or `1`-based offsets. Producers should emit this explicitly to prevent ambiguity.
* **`Dice.node(index, name; kwargs...)`**
  Associates a logical node index with a semantic name and optional extension fields.
* **`Dice.record(type::Union{Symbol, AbstractString}; kwargs...)`**
  General constructor for arbitrary key-value metadata records.
* **`Dice.metadata(obj)`**
  Serializes any arbitrary JSON-serializable object as a metadata record.

### Example Metadata Header Construction
```julia
using Dice

traj = Hybrid[
    (Int8[1, -1, 1], Float64[0.12, 0.44, 0.91]),
    (Int8[1, -1, 1], Float64[0.15, 0.43, 0.90]),
]

header = [
    Dice.comment("Simulation metadata"),
    Dice.format_record(),
    Dice.created("2026-01-26T18:44:00Z"),
    Dice.title("Trajectory for test run 42"),
    Dice.node_indexing(1),
    Dice.node(1, "CTRL_1"),
    Dice.node(3, "pass_a"),
    Dice.record(:solver; name = "my-relaxer", tolerance = 1e-8),
    Dice.metadata((git_commit = "abc123", seed = 42)),
]

Dice.save_trajectories("trajectory.dat", traj, header)
```

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
3. `state[2][free] = get_agitation(length(free))` generates a random perturbation vector and adds noise (agitation) to the continuous coordinates of unpinned nodes to help them escape local traps.
4. `Dice.propagate_pinned(model, parameters["num_steps"], state, pinned)` runs the solver dynamics for `num_steps` steps while ensuring the variables at the `pinned` indices remain unchanged.

### Example 3: Basic Model Propagation and Cut Evaluation
```julia
using Dice
using Graphs

# 1. Create a 7-vertex cycle graph
graph = Graphs.Graph(7)
add_edge!(graph, 1, 2)
add_edge!(graph, 2, 3)
add_edge!(graph, 3, 4)
add_edge!(graph, 4, 5)
add_edge!(graph, 5, 6)
add_edge!(graph, 6, 7)
add_edge!(graph, 7, 1)

# 2. Setup the model with sign-based Model II dynamics
num_steps = 100
delta_t = 4.4 / num_steps
model = Dice.Model(graph, Dice.model_2_hybrid_coupling, delta_t)

# 3. Define a specific initial hybrid state (spin_configuration, continuous_displacements)
state = (Int8[-1, 1, 1, -1, -1, 1, 1], [0.117, 0.557, -0.164, 0.122, 0.471, -0.323, 0.224])

# 4. Propagate the model for the specified number of steps
state = Dice.propagate(model, num_steps, state)

# 5. Evaluate final cut
final_cut = Dice.cut(model, state)
```
**Explanation**:
1. Constructs a standard cycle graph ($C_7$).
2. Initializes a `Dice.Model` with `Dice.model_2_hybrid_coupling` (Model II) and simulation step duration `delta_t`.
3. Defines the initial relaxed state as a `Hybrid` tuple where the first element is the discrete components (`SpinConf`) and the second is the continuous variables (`Vector{Float64}`).
4. `Dice.propagate` runs the relaxation dynamics for 100 steps, returning the final optimized hybrid state.
5. `Dice.cut` evaluates the final binary cut value of the state on the model's graph (yielding a Max-Cut size of `6`).

