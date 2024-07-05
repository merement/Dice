# A collection of functions for simulating dynamical Ising solvers
#

# NOTE:
# This is the version that phases out explicitly defined functions supporting
# the separated representation of Model II:
#     update_2!
#     step_rate_2
#     propagate_2
#     trajectories_2
#     coupling_model_2
#     cut_2

module Dice

using Arpack
using Base: Integer
using Distributions
using Graphs
using Random
using SimpleWeightedGraphs
using SparseArrays

export Hybrid, SpinConf, Model,
    dumpGraph, loadMTXAdjacency, loadMTXGraph,
    get_ER_graph,
    get_regular_graph,
    get_random_cube,
    get_random_hybrid,
    get_random_sphere,
    get_initial,
    get_random_configuration,
    sine, triangular,
    cut, get_best_cut, get_best_configuration, extract_configuration,
#    H_Ising, energy,
    number_to_conf, 
    propagate, roundup,
    test_branch, scan_for_best_configuration,
    conf_decay,
    local_search, local_twosearch,
    local_search!, local_twosearch!

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
const Hybrid = Tuple{SpinConf, FVector}

# To specify the general kind of model
const ModelGraph = Union{SimpleGraph, SimpleWeightedGraph}

const graph_types = Set([:binary, :weighted])
const isotropy_types = Set([:isotropic, :anisotropic, :dynamic])
const noise_types = Set([:noiseless, :thermal])
struct ModelKind
    graph_type::Symbol 
    anisotropy::Symbol

    Noise::Symbol

    function ModelKind()
        new(:binary, :isotropic, :noiseless)
    end
    
    function ModelKind(graph_t::Symbol)
        graph_t in graph_types || throw(ArgumentError("Invalid graph type: $graph_t"))
        new(graph_t, :isotropic, :noiseless)
    end
end

# For a compact description of simulation scenarios and controlling the
# module behavior
mutable struct Model
    # TODO: figure out the constness
    kind::ModelKind

    graph::ModelGraph
    
    coupling::Function   # dynamical coupling function
    scale::Float64       # Defines the magnitude of the timestep
                         # default = 1/max_degree(graph)

    Ks::Float64          # The magnitude of the anisotropy term
                         # default = 0
    anisotropy::Function # anisotropy function
                         # default = coupling(x, -x)

    Ns::Float64          # The magnitude of the noise term
                         # default = 0
    noise::Function      # noise function
                         # default = noiseUniform

    extended::Bool       # obsolete, phase out (inside kind)
                         # Whether there is an extension (dynamical
                         # anisotropy)
    silence::Int         # the inverse level of verbosity
                         # default = silence_default

    Model() = new(ModelKind())
    Model(x::ModelKind) = new(x)
end

# Explicit constructors
Model(graph::ModelGraph, coupling::Function) =
    begin
        M = 
        begin
            if graph isa SimpleWeightedGraph 
                Model(ModelKind(:weighted))
            else
                Model(ModelKind(:binary))
            end
        end

        M = Model()
        M.graph = graph
        M.scale = 1/Graphs.Δ(graph)
        M.coupling = coupling
#        M.method = M.coupling
        # temporary placeholder
#        M.energy = M.coupling
        # The default interpretation of anisotropy
#=         M.anisotropy = x -> M.coupling(x, -x) 
        M.Ks = 0
        M.extended = false
        M.Ns = 0
        M.noise = Dice.noiseUniform =#
        M.silence = silence_default
        return M;
    end

Model(graph::ModelGraph, coupling::Function, scale::Float64) =
    begin
        M = Model(graph, coupling)
        M.scale = scale
#=         M.Ks = 0
        M.Ns = 0
        M.noise = Dice.noiseUniform
        M.silence = silence_default
        M.extended = false
 =#        
        return M;
    end

Model(graph::ModelGraph, coupling::Function, scale::Float64, Ks::Float64) =
    begin
        M = Model(graph, coupling, scale)
        @warn "Anisotropy has no effect in the present version of Dice"
#=         M.Ks = Ks
        M.Ns = 0
        M.noise = Dice.noiseUniform
        M.extended = true         =#
        return M;
    end
        
############################################################
#
### Internal service functions
### TODO: replace by a proper logging functionality
#
############################################################

function message(model::Model, out, importance=1)
    if importance > model.silence
        println("$out ($importance)")
    end
end

const DEBUG = 1
function debug_msg(out)
    if DEBUG > 0
        println("$out (DEBUG = $DEBUG)")
    end
end

############################################################
#
### File support
#
############################################################

include("file_operations.jl")

############################################################
#
### A library of coupling functions (TODO: detach?)
#
############################################################

include("dynamical_kernels.jl")

############################################################
#
### Data processing methods
#
############################################################

"""
    roundup(V::Array{Float64})

Return `V` folded into [-2, 2].
"""
function roundup(V::Array{Float64})
    return mod.(V .+ 2, 4) .- 2
end

"""
    hybrid_to_cont(hybrid::Hybrid, r = 0.0)::Array{Float64}

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
function hybrid_to_cont(hybrid::Hybrid, r = 0.0)::Array{Float64}
    return hybrid[2] .+ hybrid[1] .+ r
end

"""
    cont_to_hybrid(V::Array{Float64,1}, 
                    r = 0.0)::Hybrid

Separate the discrete and continuous components of the given distribution
according to V = sigma + X + r, where ``sigma ∈ {-1, +1}^N``, ``X ∈ [-1, 1)^N``,
and ``-2 < r < 2`` is the rounding center.

# Arguments
    V - array to process
    r - the rounding center

OUTPUT:
    Hybrid = Tuple (sigmas, x), where
    sigmas - Int8 arrays of binary spins (-1, +1)
    xs - Float64 array of displacements [-1, 1)
"""
function cont_to_hybrid(Vinp::Array{Float64,1}, r = 0.0)::Hybrid
    V = mod.(Vinp .- r .+ 2, 4) .- 2
    sigmas = zeros(Int8, size(V))
    xs = zeros(Float64, size(V))
    # sigmas = sgn(V - r) # with sgn(0) := 1
    # xs = V - sigmas
    for i in 1:length(V)
        Vred = V[i]
        (sigmas[i], xs[i]) =
            if Vred >= 0
                (1, Vred - 1)
            else
                (-1, Vred + 1)
            end
    end
    return (sigmas, xs)
end

"""
    separate(V::Array{Float64,1}, r = 0.0)

Separate the discrete and continuous components of the given distribution
according to V = sigma + X + r, where sigma ∈ {-1, +1}^N, X ∈ [-1, 1)^N,
and -2 < r < 2 is the rounding center.

# Arguments
    V - array to process
    r - the rounding center

OUTPUT:
    sigmas - integer arrays of binary spins (-1, +1)
    xs - array of displacements [-1, 1)

    NOTE: Obsolete, see `cont_to_hybrid`
"""
function separate(Vinp::Array{Float64,1}, r = 0.0)
    V = mod.(Vinp .- r .+ 2, 4) .- 2
    sigmas = zeros(Int8, size(V))
    xs = zeros(Float64, size(V))
    # sigmas = sgn(V - r) # with sgn(0) := 1
    # xs = V - sigmas
    for i in 1:length(V)
        Vred = V[i]
        (sigmas[i], xs[i]) =
            if Vred >= 0
                (1, Vred - 1)
            else
                (-1, Vred + 1)
            end
    end
    return (sigmas, xs)
end

"""
    combine(s::Array{Int8, x::Array{Float64}, r = 0.0)

Recover the dynamic variables from the separated representation.

# Arguments
    s - the {-1, 1} array containing the discrete component
    x - the array with the continuous component
    r - the rounding center (default = 0)

NOTE: Obsolete, see `hybrid_to_cont`
"""
function combine(s::Array{Int8}, x::Array{Float64}, r = 0.0)
    return x .+ s .+ r
end


"""
    HammingD(s1::Array, s2::Array)

Evaluate the Hamming distance between binary {-1, 1}-strings `s1` and `s2` of the same length.
"""
function HammingD(s1::Array, s2::Array)
    # assert(legnth(s1) == length(s2))
    count::Int64 = 0
    for i in 1:length(s1)
        if s1[i] != s2[i]
            count += 1
        end
    end
    return count
    # return - s1 .* s2 + length(s1)    
end

"""
    EuclidD(V1::FVector, V2::FVector)

Evaluate the Euclidean distance between two distributions

# Arguments
  `V1`, `V2` - two arrays on the graph vertices

OUTPUT:
      Sum_v (V_1(v) - V_2(v))^2
"""
function EuclidD(V1::FVector, V2::FVector)
    return sum((V1 .- V2).^2)
end


############################################################
#
### Statistical methods
#
############################################################

include("statistical_methods.jl")

############################################################
#
### Cut functions
#
############################################################

"""
    cut_binary(graph::SimpleGraph, conf::SpinConf)::Int

Evaluate the cut value for the given {0, 1}-weighted graph and binary configuration

# Arguments
    graph - Graphs object
    conf - binary configuration array {-1,1}^N
OUTPUT:
    Sum_e   (1 - e1. e.2)/2
"""
function cut_binary(graph::SimpleGraph, conf::SpinConf)::Int
    if nv(graph) != length(conf)
        println("ERROR: The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match")
    end # side note: turned out to be useful

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
    Sum_e   w(e)(1 - e.1 e.2)/2
"""
function cut_weighted(graph::SimpleWeightedGraph,
                    conf::SpinConf)::Float64
    if nv(graph) != length(conf)
        println("ERROR: The configuration size $(length(conf)) and the graph size $(nv(graph)) do not match")
    end # side note: turned out to be useful

    out::Float64 = 0.0
    for edge in edges(graph)
        if conf[edge.src] * conf[edge.dst] == -1
            out += weight(edge)
        end
    end
    return out
end

"""
    cut(graph::ModelGraph, conf::SpinConf)

Main dispatch for evaluating cut. Depending on whether `graph`
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

function get_random_cut(graph, V, trials = 10)
    (becu, _, _) = get_random_rounding(graph, V, trials)
    return becu
end

function get_rate(VFull)
    # Evaluates the magnitude of variations (discrete time derivative)
    # in a 2D array VFull[time, N].
    # Returns a 1D array of magnitudes with the last element duplicated
    # to keep the same size as the number of time points in VFull

    out = [sum((VFull[:,i + 1] - VFull[:,i]).^2) for i in 1:(size(VFull)[2] - 1)]
    return sqrt.([out; out[end]])
end

# Some Model II specific functions

function cut_2(model, s::SpinConf, x::FVector)
    # This should be some kind of energy function for Model II
    # Evaluages the cut function for Model II:
    # C_2(sigma, X) = C(sigma) + \Delta C_2(sigma, X)
    # where C(\sigma) is the cut given by the binary component
    # and \Delta C_2 = \sum_{m,n} A_{m,n} s_m s_n |x_m - x_n|/2
    phix = 0
    graph = model.graph
    for edge in edges(graph)
        phix += s[edge.src]*s[edge.dst]*abs(x[edge.src] - x[edge.dst])/2
    end
    return Dice.cut(graph, s) + phix
end

function cut_2(model, distr::Hybrid)
    return cut_2(model, distr[1], distr[2])
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
    get_random_configuration(len::Int, p=0.5)

Return a Bernoulli integer sequence of length `len` and parameter `p`.
Can be used for generating random binary distributions (spin configurations).

# Arguments
    len - the length of the sequence
    p - the parameter of the Bernoulli distribution (probability to have 1)

OUTPUT:
    {-1, 1}^(len) - Int8-integer array of 1 and -1
"""
function get_random_configuration(len::Int, p=0.5)::Array{Int8}
    # A Bernoulli sequence of length len
    # Can be used for generating random binary distributions
    #    return [randspin(p) for i in 1:len]
    out = Array{Int8}(undef, len)
    return map(x ->
        if rand() < p
            1
        else
            -1
        end, out)
end

# get_random_interval
function get_initial(Nvert::Int, (vmin, vmax))
    # Generate random vector with Nvert components uniformly
    #  distributed in the (vmin, vmax) interval

    bot, top = minmax(vmin, vmax)
    mag = top - bot

    return mag .* rand(Float64, Nvert) .+ bot
end

function get_initial_2(Nvert::Int, (vmin, vmax), p=0.5)
    # Generate random configuration of length `Nvert` in the separated
    # representation with the continuous component uniformly distributed
    # in the (vmin, vmax) interval
    #
    # NOTE: OBSOLETE, replaced by get_random_hybrid

    return realign_2((get_random_configuration(Nvert, p),
                      get_initial(Nvert, (vmin, vmax))))
end

function get_random_hybrid(Nvert::Int, (vmin, vmax), p=0.5)
    # Generate random configuration of length `Nvert` in the separated
    # representation with the continuous component uniformly distributed
    # in the (vmin, vmax) interval
    #

    return realign_2((get_random_configuration(Nvert, p),
                      get_initial(Nvert, (vmin, vmax))))
end

"""
    get_random_sphere(Nvert::Int, radius)

Return vector with `Nvert` random points uniformly distributed over the
sphere of given `radius`.

# Arguments
    Nvert::Int
    radius

OUTPUT:
    Array{Float64}[1:Nvert] - points on sphere
"""
function get_random_sphere(Nvert::Int, radius)
	X = randn(Nvert)
    X ./= sqrt.(X' * X)
    return X .* radius
end

"""
    get_random_cube(Nvert::Int, side)

Return vector with `Nvert` random points uniformly distributed insde the
cube with side length given by `side`.

# Arguments
    Nvert::Int
    side

OUTPUT:
    Array{Float64}[1:Nvert] - points inside the cube
"""
function get_random_cube(Nvert::Int, side)
    return side .* rand(Float64, Nvert)
end


function randnode(nvert)
    # Returns a random number in [1, nvert]
    return rand(tuple(1:nvert...))
end

# Generator of transformations producing all strings with the fixed
# Hamming distance from the given one
#
include("hf.jl")

############################################
#
### Local search and binary transformations
#
############################################

"""
    flipconf(conf::Array, flip::Array{Int})

Change configuration `conf` according to flips in the index array `flip`
    
# Arguments
    conf - {-1, 1}^N array containing the original string
    flip - array with indices where conf should be modified

OUTPUT:
    a string at the H-distance sum(flip) from conf
"""
function flipconf(conf::Array, flip::Array{Int})
    conf[flip] .*= -1
    return conf
end

function majority_flip!(graph, conf::SpinConf, node)
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

function majority_twoflip!(graph, conf::SpinConf, cut_edge)
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

function local_search(graph, conf::SpinConf)
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
    local_search!(graph, conf)

Enforce the node majority rule in `conf` on `graph`, while changing
`conf` in-place.

This implements the 1-opt local search. It checks the nodes whether the
number of adjacent uncut edges exceeds the number of cut edges. If it does,
the node is flipped. It produces a configuration, which does not have
configurations yielding better cut within the Hamming distance one.

# Arguments
    graph - the graph object
    conf - {-1, 1}^N - the initial configuration

OUTPUT:
    count - the total number of passes
    `conf` is displaced to a locally optimal configuration
"""
function local_search!(graph, conf::SpinConf)
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

function local_twosearch(graph, conf::SpinConf)
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
    local_twosearch!(graph, conf)

Eliminate in-place pairs in configuration `conf' breaking the edge majority
rule on `graph' and return the number of passes.

The configuration is presumed to satisfy the node majority rule.

# Arguments
- `graph` - the graph object
- `conf` - {-1, 1}^N array with the initial spin configuration
"""
function local_twosearch!(graph, conf::SpinConf)
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
    number_to_conf(number ::Int, len ::Int)::Array{Int8}

Return an Int8-array of the `number`-th {-1, 1}-configuration of a model
with `len` spins.

This is `number` in the binary representation
padded with leading zeros to length
"""
function number_to_conf(number ::Int, len ::Int)::Array{Int8}
    preconf::Array{Int8} = digits(number, base=2, pad=len) |> reverse
    return 2 .* preconf .- 1
end

function conf_to_number(conf::SpinConf)
    # Convert conf as a binary number to its decimal representation

    # TODO: This function is rarely needed, hence the dumb code.
    #        I'm not even sure it works correctly
    #        It must be checked against big-endian/little-endian
    binconf = (conf .+ 1) ./ 2

    out = 0
    for i in 1:length(binconf)
        out += binconf[i] * 2^(i - 1)
    end

    return out
end

### Graph generation
#
# These are simple wrappers for functions from Graphs ensuring
# that the graphs are connected (for Erdos-Renyi, Graphs may
# return disconnected graph).

function get_connected(Nvert, prob)
    # Generate a connected Erdos-Renyi graph with `Nvert`
    # vertices and `prob` density of edges
    # More precisely. On the set of edges of a complete graph
    # K_{Nvert}, we have a (Bernoulli process) function F,
    # which takes values 1 and 0 with probabilities prob and
    # 1 - prob, respectively.
    # The output graph is a connected subgraph of K_{Nvert} with
    # only edges where F = 1 kept in the set of edges.
    #
    # NOTE Obsolete See get_ER_graph
    cnct = false
    G = Graph()
    while !cnct
        G = erdos_renyi(Nvert, prob)
        cnct = is_connected(G)
    end
    return G
end

function get_ER_graph(Nvert::Int, prob::Float64)::SimpleGraph
    # Generate a connected Erdos-Renyi graph with `Nvert`
    # vertices and `prob` density of edges
    # More precisely. On the set of edges of a complete graph
    # K_{Nvert}, we have a (Bernoulli process) function F,
    # which takes values 1 and 0 with probabilities prob and
    # 1 - prob, respectively.
    # The output graph is a connected subgraph of K_{Nvert} with
    # only edges where F = 1 kept in the set of edges.

    while true
        G = erdos_renyi(Nvert, prob)
        is_connected(G) && return G
    end
end

function get_regular_graph(Nvert, degree)::SimpleGraph
    # Generate a random connected `degree'-regular graph with `Nvert` vertices

    while true
        G = random_regular_graph(Nvert, degree)
        is_connected(G) && return G
    end
end

##########################################################
#
###   Dynamics methods
#
##########################################################

function step_rate(graph::SimpleGraph, method::Function,
                   V::FVector)::FVector
# Evaluate ΔV (continuous representation) in a single step
# Deals with {0,1}-weighted graphs, isotropic, noiseless
    out  = zeros(Float64, size(V))
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
# Evaluate ΔV (continuous representation) in a single step
# Deals with weighted graphs, isotropic, noiseless
    out  = zeros(Float64, size(V))
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

function step_rate_hybrid(graph::SimpleGraph, coupling::Function,
                     s::SpinConf, x::FVector)::FVector
    out = zeros(Float64, size(x))
    for node in vertices(graph)
        xnode = x[node]
        for neib in neighbors(graph, node)
            out[node] += s[neib]*coupling(xnode, x[neib])
        end
    end
    out .*= s
    return out
end

function step_rate_hybrid(graph::SimpleWeightedGraph, coupling::Function,
                     s::SpinConf, x::FVector)::FVector
    out = zeros(Float64, size(x))
    for node in vertices(graph)
        xnode = x[node]
        for neib in neighbors(graph, node)
            out[node] += s[neib]*coupling(xnode, x[neib])*
                graph.weights[node, neib]
        end
    end
    out .*= s
    return out
end


function trajectories(graph::ModelGraph, tmax::Int, scale::Float64,
                      method::Function, Vini)
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

    for _ = 1:(tmax - 1)
        ΔV = scale .* step_rate(graph, method, V)
        V += ΔV
        VFull = [VFull V]
    end

    return VFull
end

function trajectories(model::Model, duration, Vini)
    # Advances the graph `(duration - 1)` steps forward
    # This is the verbose version, which returns the full dynamics
    #
    # model - the model description
    # duration - how many time points to evaluate
    # V0 - the initial conditions
    #
    # OUTPUT:    #   VFull = [V(0) V(1) ... V(duration-1)]

    return trajectories(model.graph, duration, model.scale,
                        model.coupling, Vini)
end

function trajectories(graph::ModelGraph, tmax::Int, scale::Float64,
                      mtd::Function, start::Hybrid)
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps
    # Keeps the full history of progression
    S = start[1]
    X = start[2]

    SFull = start[1]
    XFull = start[2]
    
    for _ = 1:(tmax - 1)
        DX = step_rate_2(graph, mtd, S, X).*scale
        update_hybrid!(S, X, DX)
        XFull = [XFull X]
        SFull = [SFull S]        
    end
    return (SFull, XFull)
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
    for _ = 1:(tmax - 1)
        V += scale .* step_rate(graph, method, V)
    end

    return V
end

function propagate(model::Model, tmax::Int, Vini)
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
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps
    S = start[1]
    X = start[2]
    for _ = 1:(tmax - 1)
        DX = step_rate_hybrid(graph, mtd, S, X).*scale
        update_hybrid!(S, X, DX)
    end
    return (S, X)
end

function propagate(model::Model, tmax::Int, start::Hybrid)::Hybrid
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps
    return propagate(model.graph, tmax, model.scale,
                     model.coupling, start)
end


########################################################
##
###   Extended methods treating anisotropy dynamically
##
########################################################
include("dyn_anisotropy_model.jl")


##
## Functions for the separated representation
## (mostly deprecated legacies)
##

function realign_hybrid(conf::Hybrid, r = 0.0)::Hybrid
    # Changes the reference point for the separated representation by `r`
    # according to xi - r = sigma(r) + X(r)
    # INPUT & OUTPUT:
    #     conf = (sigma, X)
    V = Dice.hybrid_to_cont(conf[1], conf[2])
    return Dice.cont_to_hybrid(V, r)
end

function realign_2(conf::Hybrid, r = 0.0)
    # Changes the reference point for the separated representation by `r`
    # according to xi - r = sigma(r) + X(r)
    # INPUT & OUTPUT:
    #     conf = (sigma, X)
    #
    # NOTE: Obsolete, replaced by realign_hybrid
#    V = Dice.hybrid_to_cont(conf[1], conf[2])
    V = Dice.hybrid_to_cont(conf)
    return Dice.cont_to_hybrid(V, r)
end

function update_2!(spins::SpinConf, xs::FVector, dx::FVector)
    # Tracing version
    # Update the continuous component (xs) by dx using the wrapping rule
    # The spin part is assumed to be Int8-array
    # Return the number of flips (tracing)
    Nvert = length(spins)
    count = 0
    for i in 1:Nvert
        # we presume that |dx[i]| < 2
        Xnew = xs[i] + dx[i]
        if Xnew > 1
            xs[i] = Xnew - 2
            spins[i] *= -1
            count += 1
        elseif Xnew < -1
            xs[i] = Xnew + 2
            spins[i] *= -1
            count += 1
        else
            xs[i] = Xnew
        end
    end
    return count
end

function propagate_2(model::Model, tmax, Sstart::SpinConf, Xstart::FVector)
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps
    #
    # NOTE: Legacy version, use dispatched propagate instead
    X = Xstart
    S = Sstart
    graph = model.graph
    scale = model.scale
    mtd = model.coupling
    # Ns = model.Ns
    # noise = model.noise
    for _ = 1:(tmax - 1)
#        DX = step_rate_hybrid(graph, mtd, S, X, Ns, noise).*scale
        DX = step_rate_hybrid(graph, mtd, S, X).*scale
        update_2!(S, X, DX)
    end
    return (S, X)
end

function trajectories_2(model::Model, tmax, Sstart::SpinConf,
                        Xstart::FVector)
    # Advances the model in the initial state (Sstart, Xstart)
    # for tmax time steps
    # Keeps the full history of progression
    #
    # NOTE: Legacy version, use dispatched trajectories instead
    X = Xstart
    S = Sstart

    XFull = Xstart
    SFull = Sstart
    
    graph = model.graph
    scale = model.scale
    mtd = model.coupling
    # Ns = model.Ns
    # noise = model.noise
    for _ = 1:(tmax - 1)
#        DX = step_rate_hybrid(graph, mtd, S, X, Ns, noise).*scale
        DX = step_rate_hybrid(graph, mtd, S, X).*scale
        update_2!(S, X, DX)
        XFull = [XFull X]
        SFull = [SFull S]        
    end
    return (SFull, XFull)
end    

############################################################
#
### Simulation
#
############################################################
include("simulations.jl")

end # end of module Dice
