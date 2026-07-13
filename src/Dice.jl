# A collection of functions for simulating dynamical Ising solvers

module Dice

using Arpack
using Base: Integer
using Distributions
using Graphs
using Random
using SimpleWeightedGraphs
using SparseArrays

export
    # Main types
    Hybrid,
    SpinConf,
    Model,
    # File operations
    save_mtx_graph,
    load_mtx_adjacency,
    load_mtx_graph,
    # Generating graphs
    get_ER_graph,
    get_regular_graph,
    # Generating random states
    get_random_cube,
    get_random_sphere,
    get_random_configuration,
    get_random_hybrid,
    # Distance metrics
    hamming_distance,
    euclid_distance,
    # Cut and rounding
    cut,
    get_best_cut,
    get_best_configuration,
    extract_configuration,
    roundup,
    number_to_conf,
    # Progression methods
    propagate,
    trajectories,
    propagate_pinned,
    trajectories_pinned,
    # local search
    local_search,
    local_search!,
    local_twosearch,
    local_twosearch!

# specialized search methods like test_branch
# are not exported but are available in simulations.jl

# (reserved for future implementations)
"Define the default level of logging verbosity"
const silence_default = 3

"The smallest non-zero weight magnitude"
const weight_eps = 1e-5

########################################
##
## The main types describing the Model
##
########################################

# Data types for dynamical variables
const FVector = Vector{Float64}
const SpinConf = Vector{Int8}
# It's not a vector of RelaxedSpin = Tuple{Spin, Interval}
# because of the performance concerns
const Hybrid = Tuple{SpinConf,FVector}

# To specify the general kind of model (performance concerns?)
const ModelGraph = Union{SimpleGraph,SimpleWeightedGraph}

mutable struct Model
    graph    ::ModelGraph

    "Dynamical coupling function"
    coupling ::Function
    "The magnitude of the timestep"
    scale    ::Float64

    "the inverse level of verbosity"
    silence::Int           # default = silence_default

    Model() = new()
end

# Explicit constructors
Model(graph::ModelGraph, coupling::Function) =
    begin
        M = Model()
        M.graph = graph
        M.scale = 1 / Graphs.Δ(graph)
        M.coupling = coupling

        M.silence = silence_default
        return M
    end

Model(graph::ModelGraph, coupling::Function, scale::Float64) =
    begin
        M = Model(graph, coupling)
        M.scale = scale

        return M
    end


### Internal service functions

function message(model::Model, out, importance=1)
    if importance > model.silence
        if importance == 1
            @debug out
        elseif importance == 2
            @info out
        elseif importance == 3
            @warn out
        else
            @error out
        end
    end
end

const DEBUG = 1
function debug_msg(out)
    @debug out
end


### Data processing methods

"""
    roundup(V ::Vector{Float64}, r::Float64 = 0.0) ::Vector{Float64}

Return `V` displaced by `r` folded into `[-2, 2]` interval.
"""
function roundup(V ::Vector{Float64}, r::Float64 = 0.0) ::Vector{Float64}
    return mod.(V .+ r .+ 2, 4) .- 2
end

"""
    hybrid_to_cont(hybrid::Hybrid, r::Float64 = 0.0)::Array{Float64}

Convert the `hybrid` (σ, X) representation to the
continuous (ξ) representation with rounding center at `r`.

# Arguments
- `hybrid::Hybrid` = (s, x), where
 `s` - the {-1, 1}-array containing the discrete component
 `x` - the array with the continuous component

- `r` - the rounding center (default = 0)

# Output
    Array{Float64}(length(x)) with ξ = σ + X + r
"""
function hybrid_to_cont(hybrid::Hybrid, r::Float64=0.0)::Array{Float64}
    return hybrid[2] .+ hybrid[1] .+ r
end

"""
    hybrid_to_cont(hybrid::Hybrid, k::Vector{Int})::Vector{Float64}

Convert the `hybrid` (σ, X) representation to the
continuous (ξ) representation with rounding center at `4 k`, as
prescribed by the definition of the relaxed spin.

# Arguments
- `hybrid::Hybrid` = (s, x), where
 `s` - the {-1, 1}-array containing the discrete component
 `x` - the array with the continuous component

- `k` - array of the period displacements for each node

# Output
    Vector{Float64}(length(x)) with ξ = σ + X + 4k
"""
function hybrid_to_cont(hybrid::Hybrid, k::Vector{Int})::Vector{Float64}
    return hybrid[2] .+ hybrid[1] .+ 4 .* k
end

"""
    hybrid_to_cont_series(hybrid_series::Vector{Hybrid}) ::Vector{FVector}

Convert
"""
function hybrid_to_cont_series(hybrid_series::Vector{Hybrid})::Vector{FVector}
    # ξ = σ + X + k (node dependent)
    k_vec = zeros(Int, length(hybrid_series[1][1]))
    xi_series::Vector{FVector} = []
    push!(xi_series, hybrid_to_cont(hybrid_series[1], k_vec))

    for i in 2:length(hybrid_series)
        cur = hybrid_series[i]
        prev = hybrid_series[i-1]
        flips = findall(x -> x != 0, cur[1] - prev[1])
        for ind in flips
            if cur[2][ind] > prev[2][ind] && cur[1][ind] > prev[1][ind]
                k_vec[ind] -= 1
            elseif cur[2][ind] < prev[2][ind] && cur[1][ind] < prev[1][ind]
                k_vec[ind] += 1
            end
        end
        push!(xi_series, hybrid_to_cont(cur, k_vec))
    end
    return xi_series
end

"""
    cont_to_hybrid(V::Vector{Float64},
                    r::Float64 = 0.0)::Hybrid

Separate the discrete and continuous components of the given distribution
according to V = sigma + X + r, where ``sigma ∈ {-1, +1}^N``, ``X ∈ [-1, 1)^N``,
and ``-2 < r < 2`` is the rounding center.

# Arguments
    V - array with continuous values of dynamical variables to process
    r - the rounding center

OUTPUT:
    Hybrid = Tuple (sigmas, xs), where
           sigmas - Int8 arrays of binary spins (-1, +1)
           xs - Float64 array of displacements [-1, 1)
"""
function cont_to_hybrid(Vinp::Vector{Float64}, r::Float64=0.0)::Hybrid
    V = mod.(Vinp .- r .+ 2, 4) .- 2
    sigmas = zeros(Int8, size(V))
    xs = zeros(Float64, size(V))
    # sigmas = sgn(V - r) # with sgn(0) := 1
    # xs = V - sigmas
    for (i, Vi) in enumerate(V)
        (sigmas[i], xs[i]) = Vi >= 0 ? (1, Vi - 1) : (-1, Vi + 1)
    end
    return (sigmas, xs)
end

"""
    hamming_distance(s1::AbstractVector, s2::AbstractVector)::Int

Evaluate the Hamming distance between vectors `s1` and `s2` of the same
length. The distance is defined as the number of component-wise differences
in `s1` and `s2`.
"""
function hamming_distance(s1::AbstractVector, s2::AbstractVector)::Int
    if length(s1) != length(s2)
        throw(DimensionMismatch("Vectors lengths don't match ($(length(s1)), $(length(s2)))"))
    end
    count = 0
    @inbounds for i in eachindex(s1, s2)
        if s1[i] != s2[i]
            count += 1
        end
    end
    return count
end

"""
    euclid_distance(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real})::Float64

Evaluate the Euclidean distance between two distributions

# Arguments
  `V1`, `V2` - two arrays on the graph vertices

OUTPUT:
      Sum_v (V_1(v) - V_2(v))^2
"""
function euclid_distance(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real})::Float64
    if length(V1) != length(V2)
        throw(DimensionMismatch("Vectors lengths don't match ($(length(V1)), $(length(V2)))"))
    end
    return sum((V1[i] - V2[i])^2 for i in eachindex(V1, V2))
end


### Cut functions

"""
    cut_binary(graph::SimpleGraph, conf::SpinConf)::Int

Evaluate the integer cut value for the given {0, 1}-weighted graph and binary
configuration

# Arguments
    graph - Graphs object
    conf - binary configuration array {-1,1}^N
OUTPUT:
    return integer Σ_e (1 - e1. e.2)/2
"""
function cut_binary(graph::SimpleGraph, conf::SpinConf)::Int
    if nv(graph) != length(conf)
        throw(DimensionMismatch("The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match"))
    end # NOTE: turned out to be useful

    out::Int = 0
    for edge in edges(graph)
        if conf[edge.src] * conf[edge.dst] == -1
            out += 1
        end
    end
    return out
end

"""
    cut_weighted(graph::SimpleWeightedGraph,
                conf::SpinConf)::Float64

Evaluate the cut value for the given weighted `graph` and binary configuration `conf`

# Arguments
    graph - Graphs object
    conf - binary configuration array {-1,1}^N
OUTPUT:
    return float Σ_e w(e)(1 - e.1 e.2)/2
"""
function cut_weighted(graph::SimpleWeightedGraph,
    conf::SpinConf)::Float64
    if nv(graph) != length(conf)
        throw(DimensionMismatch("The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match"))
    end # side note: turned out to be useful

    out::Float64 = 0.0
    for edge in edges(graph)
        if conf[edge.src] * conf[edge.dst] == -1
            out += get_weight(graph, edge.src, edge.dst)
        end
    end
    return out
end

"""
    cut(graph::ModelGraph, conf::SpinConf)

The main dispatch for evaluating cut. Depending on whether `graph`
is binary or weighted, the respective cut function is called.

# Arguments
    graph - Graphs object (`ModelGraph`)
    conf - binary configuration array {-1,1}^N

# Output
    ∑_e w(e)(1 - e.1 e.2)/2
"""
function cut(graph::ModelGraph, conf::SpinConf)
    if graph isa SimpleWeightedGraph
        return cut_weighted(graph, conf)
    else
        return cut_binary(graph, conf)
    end
end

function cut(graph::ModelGraph, state::Hybrid)
    return cut(graph, state[1])
end

function cut(model::Model, state::Hybrid)
    return cut(model.graph, state[1])
end

## (TEMP: relaxed cut methods for the V-2 model)

function relaxed_cut_correction_binary_v2(graph::SimpleGraph,
                                          state::Hybrid) ::Float64
    conf = state[1]
    if nv(graph) != length(conf)
        throw(DimensionMismatch("The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match"))
    end # side note: turned out to be useful

    corr ::Float64 = 0.0
    for edge in edges(graph)
        corr += conf[edge.src] * conf[edge.dst] *
            abs(state[2][edge.src] - state[2][edge.dst])
    end
    return corr/2
end

function relaxed_cut_correction_weighted_v2(graph::SimpleWeightedGraph,
                                          state::Hybrid) ::Float64
    conf = state[1]
    if nv(graph) != length(conf)
        throw(DimensionMismatch("The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match"))
    end # side note: turned out to be useful

    corr ::Float64 = 0.0
    for edge in edges(graph)
        u, v = edge.src, edge.dst
        corr += conf[u] * conf[v] *
            get_weight(graph, u, v) * abs(state[2][u] - state[2][v])
    end
    return corr/2
end

function relaxed_cut_correction_v2(graph::ModelGraph, state::Hybrid) ::Float64
    if graph isa SimpleWeightedGraph
        return relaxed_cut_correction_weighted_v2(graph, state)
    else
        return relaxed_cut_correction_binary_v2(graph, state)
    end
end

function relaxed_cut_binary_v2(graph::SimpleGraph, state::Hybrid) ::Float64
    conf = state[1]
    if nv(graph) != length(conf)
        throw(DimensionMismatch("The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match"))
    end # side note: turned out to be useful

    cv = cut_binary(graph, state[1])
    corr = relaxed_cut_correction_binary_v2(graph, state)
    return cv + corr
end

function relaxed_cut_weighted_v2(graph::SimpleWeightedGraph,
                                 state::Hybrid) ::Float64
    conf = state[1]
    if nv(graph) != length(conf)
        throw(DimensionMismatch("The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match"))
    end # side note: turned out to be useful

    cv = cut_weighted(graph, state[1])
    corr = relaxed_cut_correction_weighted_v2(graph, state)
    return cv + corr
end

"""
    relaxed_cut_v2(graph::ModelGraph, state::Hybrid)

Evaluate the relaxed cut for 'graph' driven by the V-2 model in 'state'.
This is a dispatch making a decision based whether the passed graph is
binary or weighted

# Arguments
    graph - Graphs object (`ModelGraph`)
    conf - binary configuration array {-1,1}^N

# Output
    ∑_e w(e)(1 - e.1 e.2)/2 + 1/2 ∑_e w(e) e.1 e.2 |X.1 - X.2|
"""
function relaxed_cut_v2(graph::ModelGraph, state::Hybrid) ::Float64
    if graph isa SimpleWeightedGraph
        return relaxed_cut_weighted_v2(graph, state)
    else
        return relaxed_cut_binary_v2(graph, state)
    end
end

"""
    relaxed_cut_v2(model::Model, state::Hybrid)

Evaluate the relaxed cut for 'model' driven by the V-2 model in 'state'.
"""
function relaxed_cut_v2(model::Model, state::Hybrid) ::Float64
    return relaxed_cut_v2(model.graph, state)
end


############################################
#
### Roundings
#
############################################

include("rounding_methods.jl")

####################################################
#
###          Preparing states
#
####################################################

# Patchy Bernoulli generator
function randspin(p=0.5)::Int8
    s = rand()
    out::Int8 = s < p ? 1 : -1
    return out
end

"""
    get_random_configuration(len::Int, p=0.5) ::SpinConf

Return a Bernoulli integer sequence of length `len` and parameter `p`.
Can be used for generating random binary distributions (spin configurations).

# Arguments
    len - the length of the sequence
    p - the parameter of the Bernoulli distribution (probability to have 1)

OUTPUT:
    {-1, 1}^(len) - Int8-integer array of 1 and -1
"""
function get_random_configuration(len::Int, p=0.5) ::SpinConf
    # A Bernoulli sequence of length len
    # Can be used for generating random binary distributions
    #    return [randspin(p) for i in 1:len]
    out::SpinConf = zeros(len)
    #out = Array{Int8}(undef, len)
    return map(x ->
            if rand() < p
                1
            else
                -1
            end, out)
end

# get_random_interval

"""
    get_random_hybrid(Nvert::Int, (vmin, vmax)::Tuple{Float64, Float64}, p=0.5) ::Hybrid

Generate random hybrid configuration of length `Nvert` with the
continuous component uniformly distributed in the interval `(vmin, vmax)`
(the uniform distribution in the cube ``[vmin, vmax]^N``) and the discrete
component being a Bernoulli process with parameter `p` (the default value
`p = 0.5` means the uniform distribution in ``{-1, 1}^N``).
"""
function get_random_hybrid(Nvert::Int,
                           (vmin, vmax)::Tuple{Float64, Float64},
                           p=0.5) ::Hybrid
    vsize = vmax - vmin
    offset = (vmax + vmin)/2
    return (get_random_configuration(Nvert, p),
            get_random_cube(Nvert, size) .+ offset)
end

"""
    get_random_hybrid(Nvert::Int, vsize::Float64, p=0.5) ::Hybrid

Generate random hybrid configuration of length `Nvert` with the continuous
components uniformly distributed in the interval `(-size/2, size/2)` as in
`get_random_cube`(@ref) (the uniform distribution in the cube ``[vmin,
vmax]^N``) and the discrete component being a Bernoulli process with
parameter `p` (the default value `p = 0.5` means the uniform distribution
in ``{-1, 1}^N``).

See also [`get_random_cube`](@ref)
"""
function get_random_hybrid(Nvert ::Int,
                           vsize ::Float64,
                           p ::Float64 = 0.5) ::Hybrid
    return (get_random_configuration(Nvert, p),
            get_random_cube(Nvert, vsize))
end

"""
    get_random_sphere(Nvert::Int, radius ::Float64) ::FVector

Return vector with `Nvert` random points uniformly distributed over the
sphere of given `radius`.

# Arguments
    Nvert::Int
    radius

OUTPUT:
    Array{Float64}[1:Nvert] - points on sphere
"""
function get_random_sphere(Nvert::Int, radius ::Float64) ::FVector
    X = randn(Nvert)
    X ./= sqrt.(X' * X)
    return X .* radius
end

"""
    get_random_cube(Nvert::Int, side::Float64) ::FVector

Return vector with `Nvert` random points uniformly distributed insde the
cube centered at the origin with side length given by `side`.

# Output
    Vector{Float64}[1:Nvert] ∈ [-side/2, side/2]^Nvert
"""
function get_random_cube(Nvert::Int, side::Float64)::FVector
    return side .* (rand(Float64, Nvert) .- 0.5)
end

function randnode(nvert)
    # Returns a random number in [1, nvert]
    return rand(tuple(1:nvert...))
end


############################################
#
### Local search and binary transformations
#
############################################

"""
    flipconf(conf::Vector{Int8}, flip::Vector{Int})

Change configuration `conf` according to flips in the index array `flip`

# Arguments
    conf - {-1, 1}^N array containing the original string
    flip - array with indices where conf should be modified

OUTPUT:
    a string at the H-distance sum(flip) from conf
"""
function flipconf(conf::SpinConf, flip::Vector{Int})
    conf[flip] .*= -1
    return conf
end

"""
    majority_flip!(graph ::{SimpleGraph, SimpleWeightedGraph}, conf::SpinConf, node)

For `graph` in state `conf`, compare the weights of cut and uncut edges
incident to `node`. If the weight of uncut edges is smaller, flip the
`node`'s spin. Return `true` if the spin was flipped.
"""
function majority_flip!(graph::SimpleGraph, conf::SpinConf, node)
    # Flips conf[node] to be of the opposite sign to the majority of its neighbors

    flip_flag = false
    tot = 0
    for j in neighbors(graph, node)
        tot += conf[node] * conf[j]
    end
    if tot > 0
        conf[node] *= -1
        flip_flag = true
    end
    return flip_flag
end

function majority_flip!(graph::SimpleWeightedGraph, conf::SpinConf, node)
    # Flips conf[node] to be of the opposite sign to the majority of its neighbors

    flip_flag = false
    tot = 0
    for j in neighbors(graph, node)
        tot += conf[node] * conf[j] * weights(graph)[node, j]
    end
    if tot > 0
        conf[node] *= -1
        flip_flag = true
    end
    return flip_flag
end

"""
    majority_twoflip!(graph::SimpleGraph, conf::SpinConf, cut_edge) :Boolean

For graph `graph` in the initial state `conf` inspect the edges adjacent
to `cut_edge`, whose end points are silently assumed to belong to different
partitions. If the cut can be improved by flipping the spins at the ends
of `cut_edge`, flip the pair in place.

Return `true`, if the flip was performed, and `false` otherwise.
"""
function majority_twoflip!(graph::SimpleGraph, conf::SpinConf, cut_edge) :Boolean
    # Flips a cut pair if the edges adjacent to the cut edge
    # form the wrong majority
    # Preserves the node-majority
    flip_flag = false
    tot = 0
    for i in neighbors(graph, cut_edge.src)
        tot += conf[cut_edge.src] * conf[i]
    end
    for i in neighbors(graph, cut_edge.dst)
        tot += conf[cut_edge.dst] * conf[i]
    end

    if tot > -2
        conf[[cut_edge.src, cut_edge.dst]] .*= -1
        flip_flag = true
    end
    return flip_flag
end


"""
    local_search(graph::ModelGraph, conf::SpinConf) ::SpinConf

Perform 1-opt local search for maximum cut of (weighted) graph `graph`
initially in spin state `conf`. Return found 1-opt optimal configuration.

NOTE: The passed initial state is changed by the functon.
"""
function local_search(graph::ModelGraph, conf::SpinConf) ::SpinConf
    # Eliminates vertices breaking the majority rule
    # Attention, it changes conf
    # While it's ideologically off, it is useful for functional
    # constructions like cut(graph, local_search(graph, conf))
    nonstop = true
    while nonstop
        nonstop = false
        for node in vertices(graph)
            nonstop |= majority_flip!(graph, conf, node)
        end
    end
    return conf
end

"""
    local_search_pinned(graph::ModelGraph,
                             conf::SpinConf,
                             list_pinned ::Vector{Int}) ::SpinConf

Perform 1-opt local search on graph `graph` in initial state `conf`
while ensuring that nodes in `list_pinned` keep their initial value.

NOTE: The initial state is changed by the function (except for the pinned
spins, of course).
"""
function local_search_pinned(graph::ModelGraph,
                             conf::SpinConf,
                             list_pinned ::Vector{Int}) ::SpinConf
    nonstop = true
    while nonstop
        nonstop = false
        for node in vertices(graph)
            node in list_pinned && continue
            nonstop |= majority_flip!(graph, conf, node)
        end
    end
    return conf
end

function local_search_pinned(graph::ModelGraph,
                             conf::SpinConf,
                             pinned::Vector{Tuple{Int64,Int8}}) ::SpinConf
    pin_indices ::Vector{Int} = []
    for pin in pinned
        conf[pin[1]] = pin[2]
        push!(pin_indices, pin[1])
    end
    return local_search_pinned(graph, conf, pin_indices)
end



"""
    local_search!(graph::ModelGraph, conf)

Enforce the node majority rule in `conf` on `graph`, while changing
`conf` in-place.

This implements the 1-opt local search. It checks the nodes whether the
number of adjacent uncut edges exceeds the number of cut edges. If it does,
the node is flipped. It produces a configuration, which does not have
configurations yielding better cut within the Hamming distance one.

# Arguments
    graph - the graph object
    conf - {-1, 1}^N - the initial configuration

# Output
    count - the total number of passes
    `conf` is changed in-place to a locally optimal configuration
"""
function local_search!(graph::ModelGraph, conf::SpinConf)
    nonstop = true
    count = 0
    while nonstop
        count += 1
        nonstop = false
        for node in vertices(graph)
            nonstop |= majority_flip!(graph, conf, node)
        end
    end
    return count
end


"""
    local_twosearch(graph::SimpleGraph, conf::SpinConf) :SpinConf

Perform the second stage in the 2-opt local search of the maximum cut
of graph `graph` in the initial spin state given in `conf`. Ensure that
no pair of spins can be flipped with improving the cut. Return the 2-opt
optimal configuration.

NOTE: the function assumes that the first stage of the search was already
done, so that `conf` is the outcome of [local_search]. Consequenly, the
function checks only cut edges. For a generic initial spin configuration,
the outcome is not necessarily 2-opt optimal.
"""
function local_twosearch(graph::SimpleGraph, conf::SpinConf) :SpinConf
    # Eliminates pairs breaking the edge-majority rule
    # Attention, it changes conf
    nonstop = true
    while nonstop
        nonstop = false
        for link in edges(graph)
            if conf[link.src] * conf[link.dst] < 1
                nonstop |= majority_twoflip!(graph, conf, link)
            end
        end
    end
    return conf
end

"""
    local_twosearch!(graph::SimpleGraph, conf)

Eliminate in-place pairs in configuration `conf' breaking the edge majority
rule on `graph' and return the number of passes.

The configuration is presumed to satisfy the node majority rule.

# Arguments
- `graph` - the graph object
- `conf` - {-1, 1}^N array with the initial spin configuration
"""
function local_twosearch!(graph::SimpleGraph, conf::SpinConf)
    nonstop = true
    count = 0
    while nonstop
        count += 1
        nonstop = false
        for link in edges(graph)
            if conf[link.src] * conf[link.dst] < 1
                nonstop |= majority_twoflip!(graph, conf, link)
            end
        end
    end
    return count
end

"""
    number_to_conf(number ::Int, len ::Int)::SpinConf

Return an Int8-array of the `number`-th {-1, 1}-configuration of a model
with `len` spins.

This is `number` in the binary representation padded with leading zeros to
length
"""
function number_to_conf(number::Int, len::Int)::SpinConf
    preconf::SpinConf = digits(number, base=2, pad=len) |> reverse
    return 2 .* preconf .- 1
end

"""
    conf_to_number(conf::SpinConf)::Int

Evaluate the number of configuration `conf`.

The number is obtained by treating `conf` as a binary number and converting
it to its decimal representation. Inverse of [`number_to_conf!`](@ref).
"""
function conf_to_number(conf::SpinConf)::Int
    str_form = map(x -> x == 1 ? "1" : x == -1 ? "0" : "X", conf) |> join
    return parse(Int, str_form, base=2)
end

### Graph generation
#
# These are simple wrappers for functions from Graphs ensuring
# that the graphs are connected (for Erdos-Renyi, Graphs may
# return disconnected graph).

"""
    get_ER_graph(Nvert::Int, prob::Float64)::SimpleGraph

Generate a connected Erdos-Renyi graph with `Nvert`
vertices and `prob` density of edges.

On the set of edges of a complete graph ``K_{Nvert}``, we have a (Bernoulli
process) function ``F`` taking values 1 and 0 with probabilities ``p`` and
``1 - p``, respectively. The output graph is a connected subgraph of
``K_{Nvert}`` with only edges where ``F = 1`` kept in the set of edges.
"""
function get_ER_graph(Nvert::Int, prob::Float64)::SimpleGraph
    while true
        G = erdos_renyi(Nvert, prob)
        is_connected(G) && return G
    end
end


"""
    get_regular_graph(Nvert, degree)::SimpleGraph

Generate connected random `degree`-regular {0, 1}-weighted graph.
"""
function get_regular_graph(Nvert, degree)::SimpleGraph
    # Generate a random connected `degree'-regular graph with `Nvert` vertices
    while true
        G = random_regular_graph(Nvert, degree)
        is_connected(G) && return G
    end
end

###   Dynamics methods

function step_rate(graph::SimpleGraph, method::Function,
    V::FVector)::FVector
    # Evaluate Δξ (ontinuous representation) for a single step
    # NOTE: for {0,1}-weighted graphs, isotropic, noiseless
    out = zeros(Float64, size(V))
    for node in vertices(graph)
        Vnode = V[node]
        for neib in neighbors(graph, node)
            out[node] += method(Vnode, V[neib])
        end
    end
    return out
end

function step_rate(graph::SimpleWeightedGraph, method::Function,
    V::FVector)::FVector
    # Evaluate Δξ (continuous representation) for a single step
    # NOTE: for {0,1}-weighted graphs, isotropic, noiseless
    out = zeros(Float64, size(V))
    for node in vertices(graph)
        Vnode = V[node]
        for neib in neighbors(graph, node)
            out[node] += method(Vnode, V[neib]) * graph.weights[node, neib]
        end
    end
    return out
end

"""
    update_hybrid!(spins::SpinConf, xs::FVector, dx::FVector)

Update the continuous component `xs` by `dx` using the wrapping rule with
the spin component in `spins`. Assume that |`dx[i]`| < 2.
"""
function update_hybrid!(spins::SpinConf, xs::FVector, dx::FVector)
    Nvert = length(spins)
    Xnew = xs + dx
    for i in 1:Nvert
        if Xnew[i] > 1
            xs[i] = Xnew[i] - 2
            spins[i] *= -1
        elseif Xnew[i] < -1
            xs[i] = Xnew[i] + 2
            spins[i] *= -1
        else
            xs[i] = Xnew[i]
        end
    end
end

function step_rate_hybrid(graph::SimpleGraph,
                          coupling::Function,
                          s::SpinConf,
                          x::FVector)::FVector
    out = zeros(Float64, size(x))
    for node in vertices(graph)
        xnode = x[node]
        for neib in neighbors(graph, node)
            out[node] += s[neib] * coupling(xnode, x[neib])
        end
    end
    out .*= s
    return out
end

function step_rate_hybrid(graph::SimpleWeightedGraph,
                          coupling::Function,
                          s::SpinConf,
                          x::FVector)::FVector
    out = zeros(Float64, size(x))
    for node in vertices(graph)
        xnode = x[node]
        for neib in neighbors(graph, node)
            out[node] += s[neib] * coupling(xnode, x[neib]) *
                         graph.weights[node, neib]
        end
    end
    out .*= s
    return out
end

"""
    step_rate_hybrid_pinned(graph::SimpleWeightedGraph,
                            coupling::Function,
                            s::SpinConf,
                            x::FVector,
                            pinned::Vector{Int64})::FVector

Evaluate the vector of rates for a weighted `graph` with the coupling
function `coupling`, discrete spin state `s` and continuous component `x`.
This function presumes that there are pinned nodes, for which the rate
is set to zero.
"""
function step_rate_hybrid_pinned(graph::SimpleWeightedGraph,
                                 coupling::Function,
                                 s::SpinConf,
                                 x::FVector,
                                 pinned::Vector{Int64})::FVector
    out = zeros(Float64, size(x))
    for node in vertices(graph)
        node in pinned && continue

        for neib in neighbors(graph, node)
            out[node] += s[neib] * coupling(x[node], x[neib]) *
                         graph.weights[node, neib]
        end
    end
    out .*= s
    return out
end


function trajectories(graph::ModelGraph, tmax::Int, scale::Float64,
    method::Function, Vini ::Vector{Float64}) ::Vector{Vector{Float64}}
    # Advances the graph duration - 1 steps forward
    # This is the verbose version, which returns the full dynamics
    #
    # scale - parameter to tweak the dynamics
    # duration - how many time points to evaluate
    # V0 - the initial conditions
    #
    # OUTPUT:
    #   VFull = [V(0) V(1) ... V(duration-1)]

    VFull = Vini
    V = Vini

    for _ = 1:(tmax-1)
        ΔV = scale .* step_rate(graph, method, V)
        V += ΔV
        VFull = [VFull V]
    end

    return VFull
end

function trajectories(graph::ModelGraph, tmax::Int, scale::Float64,
    mtd::Function, start::Hybrid)::Vector{Hybrid}
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps
    # Keeps the full history of progression
    S = start[1]
    X = start[2]

    state_full::Vector{Hybrid} = []

    for _ = 1:(tmax-1)
        DX = step_rate_hybrid(graph, mtd, S, X) .* scale
        update_hybrid!(S, X, DX)
        push!(state_full, (copy(S), copy(X)))
    end
    return state_full
end

function trajectories(model::Model, tmax::Int, start::Hybrid)
    return trajectories(model.graph, tmax, model.scale,
        model.coupling, start)
end


function propagate(graph::ModelGraph, tmax::Int, scale::Float64,
    method::Function, Vini::FVector)::FVector
    # Advances the graph duration - 1 steps forward
    # This is the short version, which returns only the final state vector
    #
    # scale - parameter to tweak the dynamics
    # duration - how many time points to evaluate
    # Vini - the initial conditions
    #
    # OUTPUT:
    #   [V[1] .. V[nv(graph)] at t = duration - 1
    #
    # NOTE: the order of parameters changed in 0.2.0

    V = Vini
    for _ = 1:(tmax-1)
        V += scale .* step_rate(graph, method, V)
    end

    return V
end

function propagate(model::Model, tmax::Int, Vini::FVector)
    # Advances the model::Model duration - 1 steps forward
    # Returns only the final state vector
    #
    # model - the model description
    # duration - how many time points to evaluate
    # Vini - the initial conditions
    #
    # OUTPUT:
    #   [V[1] .. V[nv(graph)] at t = duration - 1
    return propagate(model.graph, tmax, model.scale,
        model.coupling, Vini)
end

function propagate(graph::ModelGraph, tmax::Int, scale::Float64,
    mtd::Function, start::Hybrid)::Hybrid
    # Advance the model in the initial hybrid state (Sstart, Xstart)
    # for tmax time steps
    S = start[1]
    X = start[2]
    for _ = 1:(tmax-1)
        DX = step_rate_hybrid(graph, mtd, S, X) .* scale
        update_hybrid!(S, X, DX)
    end
    return (S, X)
end

function propagate(model::Model, tmax::Int, start::Hybrid)::Hybrid
    # Advance `model` in the hybrid initial state `start = (Sstart, Xstart)`
    # for `tmax` time steps
    return propagate(model.graph, tmax, model.scale,
        model.coupling, start)
end

### (TEMP) Methods propagating the system with pinned spins

"""
    propagate_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                          mtd::Function, start::Hybrid,
                          pinned::Vector{Tuple{Int64, Int8}})::Hybrid

Advance the network set on `graph` governed by the dynamical kernel `mtd`
in the initial state `start` for `tmax` steps each of duration `scale` with
spins listed in `pinned` together with the presumed values fixed.

NOTE: This version ensures that the discrete component of the pinned
relaxed spins is set to the value specified in `pinned`. The continuous
component of the pinned nodes is set to 0.0. For a version of the method that
simply preserves the initial values of the pinned relaxed spins, see the method
with `pinned` being `Vector{Int64}`.
"""
function propagate_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
    mtd::Function, start::Hybrid,
    pinned::Vector{Tuple{Int64,Int8}})::Hybrid

    S = start[1]
    X = start[2]
    # ensure that pins are set correctly for the interface uniformity
    for pin in pinned
        S[pin[1]] = pin[2]
        X[pin[1]] = 0.0
    end
    pin_pos = [pin[1] for pin in pinned]
    for _ = 1:(tmax-1)
        DX = step_rate_hybrid_pinned(graph, mtd, S, X, pin_pos) .* scale
        update_hybrid!(S, X, DX)
        # This code is kept here to remind that the pinned version of
        # step_rate explicitly avoids updating pinned nodes
        # for pin in pinned
        #     S[pin[1]] = pin[2]
        #     X[pin[1]] = 0.0
        # end
    end
    return (S, X)
end

"""
    propagate_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                          mtd::Function, start::Hybrid,
                          pinned::Vector{Int64})::Hybrid

Advance the network set on `graph` governed by the dynamical kernel `mtd`
in the initial state `start` for `tmax` steps each of duration `scale` with
spins listed in `pinned` fixed.

NOTE: This version preserves the values of pinned relaxed spins (both the
binary and continuous components) as supplied in `state`.
"""
function propagate_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                          mtd::Function, start::Hybrid,
                          pinned::Vector{Int64})::Hybrid

    S = start[1]
    X = start[2]
    # # TEMP: ensure that pinned pins are preserved
    # pinned_spins = [(S[pin], X[pin]) for pin in pinned]

    for _ = 1:(tmax-1)
        DX = step_rate_hybrid_pinned(graph, mtd, S, X, pinned) .* scale
        update_hybrid!(S, X, DX)
    end
    # # TEMP: verify the pinned nodes
    # for (ind, pin) in enumerate(pinned)
    #     if S[pin] != pinned_spins[ind][1] ||
    #         abs(X[pin] - pinned_spins[ind][2]) > weight_eps
    #         error("Pinned spin #$i (index $pin) is not preserved")
    #     end
    # end
    return (S, X)
end

function propagate_pinned(model::Model, tmax::Int, start::Hybrid,
                          pinned::Vector)::Hybrid
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps with pinned spins in pinned
    # Should work for both pinned versions
    return propagate_pinned(model.graph, tmax, model.scale,
        model.coupling, start, pinned)
end

# function propagate_pinned(model::Model, tmax::Int, start::Hybrid,
#                           pinned::Vector{Tuple{Int64,Int8}})::Hybrid
#     # Advances the model in the initial state (Sstart, Xstart)
#     # for tmax time steps with pinned spins in pinned
#     return propagate_pinned(model.graph, tmax, model.scale,
#         model.coupling, start, pinned)
# end

"""
    trajectories_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                        mtd::Function, start::Hybrid,
                        pinned::Vector{Tuple{Int64, Int8}}) ::Vector{Hybrid}

Advance the network set on `graph` governed by the dynamical kernel `mtd`
in the initial relaxed spin state `start` for `tmax` steps each of duration
`scale` with spins listed in `pinned` fixed. Store all intermediates states
of the relaxed spin vector at all instances, including the initial state, and
return the accumulated array of relaxed spin states.

NOTE: This version ensures that the discrete component of the pinned
relaxed spins is set to the value specified in `pinned`. The continuous
component of the pinned nodes is set to 0.0. For a version of the method that
simply preserves the initial values of the pinned relaxed spins, see the method
with `pinned` being `Vector{Int64}`.
"""
function trajectories_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                             mtd::Function, start::Hybrid,
                             pinned::Vector{Tuple{Int64,Int8}})::Vector{Hybrid}
    S = start[1]
    X = start[2]
    for pin in pinned
        S[pin[1]] = pin[2]
        X[pin[1]] = 0.0
    end
    pin_pos = [pin[1] for pin in pinned]

    sx::Vector{Hybrid} = []
    push!(sx, (copy(S), copy(X)))

    for _ = 1:(tmax-1)
        DX = step_rate_hybrid_pinned(graph, mtd, S, X, pin_pos) .* scale
        update_hybrid!(S, X, DX)
        push!(sx, (copy(S), copy(X)))
    end
    return sx
end

"""
    trajectories_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                        mtd::Function, start::Hybrid,
                        pinned::Vector{Int64}) ::Vector{Hybrid}

Advance the network set on `graph` governed by the dynamical kernel `mtd`
in the initial relaxed spin state `start` for `tmax` steps each of duration
`scale` with spins listed in `pinned` fixed. Store all intermediates states
of the relaxed spin vector at all instances, including the initial state, and
return the accumulated array of relaxed spin states.

NOTE: This version preserves the initial value of pinned relaxed spins.
"""
function trajectories_pinned(graph::ModelGraph, tmax::Int, scale::Float64,
                             mtd::Function, start::Hybrid,
                             pinned::Vector{Int64})::Vector{Hybrid}
    S = start[1]
    X = start[2]
    # # TEMP: ensure that pinned pins are preserved
    # pinned_spins = [(S[pin], X[pin]) for pin in pinned]

    sx::Vector{Hybrid} = []
    push!(sx, (copy(S), copy(X)))

    for _ = 1:(tmax-1)
        DX = step_rate_hybrid_pinned(graph, mtd, S, X, pinned) .* scale
        update_hybrid!(S, X, DX)
        push!(sx, (copy(S), copy(X)))
    end
    # # TEMP: verify the pinned nodes
    # for (ind, pin) in enumerate(pinned)
    #     if S[pin] != pinned_spins[ind][1] ||
    #         abs(X[pin] - pinned_spins[ind][2]) > weight_eps
    #         error("Pinned spin #$i (index $pin) is not preserved")
    #     end
    # end

    return sx
end

function trajectories_pinned(model::Model, tmax::Int, start::Hybrid,
                             pinned::Vector)::Vector{Hybrid}
    return trajectories_pinned(model.graph, tmax, model.scale,
                               model.coupling, start, pinned)
end

# function trajectories_pinned(model::Model, tmax::Int, start::Hybrid,
#                              pinned::Vector{Tuple{Int64,Int8}})
#     return trajectories_pinned(model.graph, tmax, model.scale,
#                                model.coupling, start, pinned)
# end


### Functions for the separated (relaxed spin) representation

function realign_hybrid(conf::Hybrid, r=0.0)::Hybrid
    # Changes the reference point for the separated representation by `r`
    # according to xi - r = sigma(r) + X(r)
    # INPUT & OUTPUT:
    #     conf = (sigma, X)
    # Rotates the circle by r, and increases all continous components by -r
    V = Dice.hybrid_to_cont(conf, 0.0)
    return Dice.cont_to_hybrid(V, r)
end

### Simulation

include("simulations.jl")

### File support

include("file_operations.jl")

### Statistical methods

include("statistical_methods.jl")


### A library of coupling functions

include("dynamical_kernels.jl")

### Assorted deprecated methods

include("deprecated.jl")

end # end of module Dice
