# Dynamical Ising solver
#
# Here, the dynamical variables are regarded from the perspective of the
# feasible configurations. Therefore, the Ising states are {-1, 1}^N,
# and the period of the main functions is [-2, 2).

# The main functions are F(v) and f = F'(v)
#
# F(v) is defined by
#
# C_G(v) = \sum_{(m,n) \in E} F(v_m - v_n), when v_m ∈ {-1, 1}
#
# So that F(0) = 0 (the same partition), F(±2) = -1. Hence, the period above.

# TODO:
#   1. Make energy methods for correct energy evaluations
#      Support relaxed objective functions 
#   2. Enable LUT-defined methods
#   3. Add the weight support
#   4. Implement logging
#   5. Fix the chain of types
#   6. Admit differenttial equations 

module Dice

using Base:Integer
using Arpack
using LightGraphs
using SparseArrays

export Model,
    get_connected,
    get_initial,
    sine, triangular,
    cut, get_best_cut, get_best_configuration, extract_configuration,
    number_to_conf, 
    propagate, roundup,
    test_branch, scan_vicinity, scan_for_best_configuration,
    conf_decay,
    local_search, local_twosearch

# The main type describing the Model
# For a compact description of simulation scenarios and controlling the
# module behavior
const silence_default = 3
mutable struct Model
    graph::SimpleGraph
    method::Function
    # method_energy::Function
    scale::Float64  # Defines the magnitude of the exchange terms
    Ks::Float64     # The magnitude of the anisotropy term in scale's (not used)
    silence::Integer # the inverse level of verbosity
    Model(graph, method, scale) = new(graph, method, scale, 0, silence_default)
    Model(graph, method, scale, Ks) = new(graph, method, scale, Ks, silence_default)
    Model(graph, method, scale, Ks, silence) = new(graph, method, scale, Ks, silence)
end

############################################################
#
### Internal service functions
#
############################################################

function message(model::Model, out, importance=1)
    # must be replaced by a proper logging functionality
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

"""
Read graph in the "reduced" Matrixmarket format.

Usage:
    G = loadDumpedGraph(filename::String)::SimpleGraph

INPUT:
    filename::String name of the file

OUTPUT:
    SimpleGraph object

The reduced format:
    Header: lines starting with %
    |V| |E|
    u v w_{u,v}
"""
function loadDumpedGraph(filename::String)::SimpleGraph
    # TODO simple validation

    mtxFile = open(filename, "r")
    tokens = split(lowercase(readline(mtxFile)))
    # Skip commenting header
    while tokens[1][1] == '%'
        global tokens
        tokens = split(lowercase(readline(mtxFile)))
    end

    Nodes, Edges = map((x) -> parse(Int64, x), tokens)
    G = LightGraphs.SimpleGraph(Nodes)

    for line in readlines(mtxFile)
        global tokens
        global countEdges
        tokens = split(lowercase(line))
        # we drop the weight for now
        u, v = map((x) -> parse(Int64, x), tokens[1:2])
        add_edge!(G, u, v)
    end
    close(mtxFile)

    return G
end

"""
    dumpGraph(Graph::SimpleGraph, filename::String)

Writes the sparce adjacency matrix of `Graph` to file `filename`.

NOTE: the current version writes the (0,1) matrix only.

The adjacency matrix is written in the reduced `MatrixMarket` format
format. The first line contains two numbers: `|V| |E|`
Next lines contain three numbers: `u v A_{u, v}`
where `u` and `v` are the nodes numbers and `A_{u, v}` is the edge weight.
"""
function dumpGraph(G, filename)
    weight = 1
    open(filename, "w") do out
        println(out, LightGraphs.nv(G), " ", LightGraphs.ne(G))
        for edge in LightGraphs.edges(G)
            println(out, edge.src, " ", edge.dst, " ", weight)
        end
    end
end

############################################################
#
# Common methods for evaluating the change rate and the Hamiltonian
#
############################################################

function none_method(v)
    return 0 * v  # for the type stability thing
end

function none_method(v1, v2)
    return 0 * v1  # for the type stability thing
end

const Pi2 = pi / 2
const Pi4 = pi / 4

"""
    cosine(v)

The coupling energy inducing the XY-model (rank-2 relaxation).
"""
function cosine(v)
    return cos.(Pi2 .* v)
end

function sine(v)
    return sin.(Pi2 .* v)
end

function sine(v1, v2)
    return sine(v1 - v2)
end

function linear(v1, v2)
    return v1 - v2
end

function dtri(v)
    # the derivative of the triangular function
    p = Int.(floor.(v ./ 2 .+ 1 / 2))
    parity = rem.(abs.(p), 2)
    return -2 .* parity .+ 1
end

function triangular(v)
    # local p::Integer
    # local parity::Integer
    # local envelope::Integer
    p = Int.(floor.(v ./ 2 .+ 1 / 2))
    return (v .- p .* 2) .* dtri(v)
end

function triangular(v1, v2)
    return triangular(v1 - v2)
end

const pwDELTA = 0.1
# const pwPERIOD = 4.0

function piecewise(v)
    Delta = pwDELTA
    vbar = mod.(v .+ 2, 4) .- 2
    s = sign.(vbar)

    ind = sign.(s .* vbar .- 1)
    out = Delta .* ind .+ 0.5

    return s .* out
end

function piecewise(v1, v2)
    return piecewise(v1 - v2)
end

function square(v)
    return dtri(v .+ 1)
end

function square(v1, v2)
    return square(v1 - v2)
end

function bilinear(v1, v2)
    return -dtri(v1) .* triangular(v2)
end

function bisine(v1, v2)
    return -Pi2 .* cos(Pi2 .* v1) .* triangular(Pi2 .* v2)
end

function squarishk(v, k=10)
    return tanh.(k .* triangular(v))
end

function squarish(v1, v2, k=10)
    return squarishk(v1 - v2, k)
end

const stW = 0.55 * 2

function skew_triangular(v)
    c1 = 1 / (stW * (stW - 2))

    vbar = mod.(v .+ 2, 4) .- 2
    s = sign.(vbar)
    svbar = s .* vbar
    ind = sign.(svbar .- stW)

    mid = (svbar .* (stW - 1) .- stW) .* c1
    Delta = (svbar .- stW) .* c1

    out = Delta .* ind .+ mid
    return s .* out
end

function skew_triangular(v1, v2)
    return skew_triangular(v1 - v2)
end

############################################################
#
### Analysis methods
#
############################################################

function roundup(V)
    # returns V folded into the interval [-2, 2]

    return mod.(V .+ 2, 4) .- 2
end

function HammingD(s1, s2)
    # Evaluates the Hamming distance between binary strings s1 and s2

    count = 0
    for i in 1:length(s1)
        if s1[i] != s2[i]
            count += 1
        end
    end
    return count
end

function EuclidD(V1, V2)
    # Evaluates the Euclidean distance between two distributions
    #
    # INPUT:
    #   V1, V2 - two arrays on the graph vertices
    #
    # OUTPUT:
    #   Σ_v (V_1(v) - V_2(v))^2

    return sum((V1 .- V2).^2)
end

function av_dispersion(graph, V)
    # Evaluates average difference between the dynamical variables
    # at adjacent nodes
    #
    # INPUT:
    #   graph
    #   V - the array of dynamical variables
    #
    # OUTPUT:
    #   Σ_e |V(a) - V(b)|/|graph|

    return sum([abs(V[edge.src] - V[edge.dst]) for edge in edges(graph)]) / nv(graph)
end

function c_variance(V, intervals)
    # Calculates first two momenta of the pieces of data that fall within given M intervals
    #
    # INPUT:
    #   V - array of data
    #   intervals - M x 2 array of intervals boundaries
    #
    # OUTPUT:
    #   M x 3 array with the number of points, mean and variance of data inside the respective intervals

    out = zeros(size(intervals)[1], size(intervals)[2] + 1)

    for i in 1:(size(out)[1])
        ari = filter(t -> intervals[i,1] < t < intervals[i,2], V)
        ni = length(ari)
        av = sum(ari) / ni
        var = sum((ari .- av).^2) / ni
        out[i, 1] = ni
        out[i, 2] = av
        out[i, 3] = sqrt(var)
    end

    return out
end

# function H_Ising(graph, conf)
#     # Evaluates the Ising energy for the given graph and configuration
#     #
#     # INPUT:
#     #   graph - LightGraphs object
#     #   conf - configuration array with elemnts \pm 1
#     #
#     # OUTPUT:
#     #   conf * A * conf /2 energy

#     en = 0
#     for edge in edges(graph)
#         en += conf[edge.src]*conf[edge.dst]
#     end
#     return en
# end

"""
    cut(graph, conf)

Evaluate the (weighted) cut value for the given graph and binary configuration

INPUT:
    graph - LightGraphs object
    conf - binary configuration array with elemnts \pm 1
OUTPUT:
    sum_e (1 - e1. e.2)/2
"""
function cut(graph, conf)

    if nv(graph) != length(conf)
        num = conf_to_number(conf)
        println("ERROR: The configuration $num and the vertex set have different size")
    end

    out = 0
    for edge in edges(graph)
        out += (1 - conf[edge.src] * conf[edge.dst]) / 2
    end
    return out
end

function get_rate(VFull)
    # Evaluates the magnitude of variations (discrete time derivative)
    # in a 2D array VFull[time, N].
    # Returns a 1D array of magnitudes with the last element duplicated
    # to keep the same size as the number of time points in VFull

    out = [sum((VFull[:,i + 1] - VFull[:,i]).^2) for i in 1:(size(VFull)[2] - 1)]
    return sqrt.([out; out[end]])
end

function get_best_rounding(graph, V)
    # Finds the best cut following the CirCut algorithm
    # Returns the cut, the configuration, and the threshold (t_c)
    #
    # NOTE:
    #   The algorithm operates with the left boundary of the identifying
    #   interval. The function extract_configuration, in turn, asks for the
    #   rounding center (bad design?).
    #
    # INPUT:
    #   graph
    #   V - array with the voltage distribution assumed to be within [-2, 2]
    #
    # OUTPUT:
    #   (bestcut, bestconf, bestbnd)
    #           where
    #               bestcut - the best cut found
    #               bestconf - rounded configuration
    #               bestbnd - the position of the rounding center (t_c)

    Nvert = nv(graph)

    bestcut = -1
    bestth = -2
    bestconf = zeros(Nvert)

    d = 2 # half-width of the interval

    vvalues = sort(V)
    push!(vvalues, 2)
    threshold = -2

    start = 1
    stop = findfirst(t -> t > 0, vvalues)
    if isnothing(stop) # there's no such element
        stop = Nvert + 1
    end

    while threshold < 0
        # here, we convert the tracked left boundary to the rounding center
        conf = extract_configuration(V, threshold + 1)
        cucu = cut(graph, conf)
        if cucu > bestcut
            bestcut = cucu
            bestth = threshold
            bestconf = conf
        end
        if vvalues[start] <= vvalues[stop] - d
            threshold = vvalues[start]
            start += 1
        else
            threshold = vvalues[stop] - d
            stop += 1
        end
    end
    # bestth is converted since the center is expected
    return (bestcut, bestconf, bestth + 1)  
end

function get_best_configuration(graph, V)
    # Finds the best cut following the CirCut algorithm
    # Returns the cut and the configuration
    #
    # INPUT:
    #   graph
    #   V - array with the voltage distribution assumed to be within [-2, 2]
    #
    # OUTPUT:
    #   (bestcut, bestbnd)
    #           where
    #               bestcut - the best cut found
    #               bestconf - rounded configuration

    (becu, beco, beth) = get_best_rounding(graph, V)
    return (becu, beco)
end

function get_best_cut(graph, V)
    # Finds the best cut following the CirCut algorithm
    # Returns the cut
    #
    # INPUT:
    #   graph
    #   V - array with the voltage distribution assumed to be within [-2, 2]
    #
    # OUTPUT:
    #   bestcut - the best cut found

    (becu, beco, beth) = get_best_rounding(graph, V)
    return becu
end

####################################################
#
###          Service methods
#
####################################################

# Patchy Bernoulli generator
function randspin(p=0.5)
    s = rand()
    return s < p ? 1 : -1
end

function randvector(len::Integer, p=0.5)
    # A Bernoulli sequence of length len
    return [randspin(p) for i in 1:len]
end

function randnode(nvert)
    # Returns a random number in [1, nvert]
    return rand(tuple(1:nvert...))
end

# Generator of transformations producing all strings with the fixed
# Hamming distance from the given one
#
# USAGE: hamseq[depth](length)
# `depth` is the Hamming distance
# `length` is the bitwise length of the string
#
# OUTPUT is [depth, C]-array of flip indices,
# where C is the number of strings at the given HD
#
# NOTE: Currently `depth` is limited by 5
# TODO: make a universal version
include("hf.jl")

function flipconf(conf, flip)
    # Changes configuration conf according to flips in the index array flip
    #
    # INPUT:
    #   conf - {+1, -1}^N array containing the original string
    #   flip - array with indices where conf should be modified
    #
    # OUTPUT:
    #   a string at the H-distance length(flip) from conf
    #
    # Q: isn't this conf[flip] .*= -1?

    for ind in flip
        conf[ind] *= -1
    end
    return conf
end

function majority_flip!(graph, conf, node)
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

function majority_twoflip!(graph, conf, cut_edge)
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
        conf[cut_edge.src] *= -1
        conf[cut_edge.dst] *= -1
        flip_flag = true
    end
    return flip_flag
end

function local_search(graph, conf)
    # Eliminates vertices breaking the majority rule
    nonstop = true
    while nonstop
        nonstop = false
        for node in vertices(graph)
            nonstop |= majority_flip!(graph, conf, node)
        end
    end
    return conf
end

function local_twosearch(graph, conf)
    # Eliminates pairs breaking the edge-majority rule
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
    number_to_conf(number, length)

Return `number`-th configuration of a model with `length` spins

This is `number` in the binary representation
padded with leading zeros to length
"""
function number_to_conf(number, length)
    preconf = digits(number, base=2, pad=length) |> reverse
    return 2 .* preconf .- 1
end

function conf_to_number(conf)
    # Converts conf as a binary number to its decimal representation

    # TODO: This function is rarely needed, hence the dumb code.
    #        I'm not even sure it works correctly
    #        It must be checked against big-endian/little-endian
    binconf = (conf .+ 1) ./ 2

    out = 0
    for i = 1:length(binconf)
        out += binconf[i] * 2^(i - 1)
    end

    return out
end

function extract_configuration(V::Array, threshold) # , minimal = false
    # Binarizes V according to the threshold
    # In the modular form, the mapping looks like
    #
    # V ∈ [threshold, threshold + width] -> C_1
    # V ∈ [threshold - width, threshold] -> C_2
    #
    # OPTION: On top of this, we use the global sign inversion
    #          symmetry (disabled, untested)
    #
    # INPUT:
    #   V - data array (is presumed to be rounded and within [-2, 2])
    #   threshold - the rounding center
    #   width - (TODO: the width of the central interval)
    #
    # OUTPUT:
    #   size(V) array with elements + 1 and -1 depending on the relation of
    #           the respective V elements with threshold

    width = 1 # half-width of the rounding interval

    # if sum(abs.(V))/length(V) > 2
    #     println("Error: V value is out of bounds")
    # end

    if abs(threshold) <= 1
        inds = threshold - width .<= V .< threshold + width
    else
        return -extract_configuration(V, threshold - 2 * sign(threshold))
    end
    out = 2 .* inds .- 1

    # # if we want the outcome with smaller total displacement
    # if minimal && sum(abs.(V .+ out)) < sum(abs.(V .- out))
    #     out .*= -1
    # end
    return out
end

function get_connected(Nvert, prob)
    # Generates a connected graph with `Nvert` vertices and
    # `prob` density of edges
    # A bit more precisely. On the set of edges of a complete graph
    # K_{Nvert}, we have a (Bernoulli process) function F, which takes
    # values 1 and 0 with probabilities prob and 1 - prob, respectively.
    # The output graph is a connected subgraph of K_{Nvert} with
    # only edges where F = 1 kept in the set of edges.
    cnct = false
    G = Graph()
    while !cnct
        G = erdos_renyi(Nvert, prob)
        cnct = is_connected(G)
    end
    return G
end

function get_initial(Nvert::Integer, (vmin, vmax))
    # Generate random vector with Nvert components uniformly distributed
    # in the (vmin, vmax) interval

    bot, top = minmax(vmin, vmax)
    mag = top - bot

    return mag .* rand(Float64, Nvert) .+ bot
end


##########################################################
#
###   Dynamics methods
#
##########################################################


function step_rate(graph::SimpleGraph, methods::Array{Function}, V::Array, Ks::Float64)
    # Evaluates ΔV for a single step
    # This version takes an array of methods evaluating coupling
    # and local energies
    #
    # INPUT:
    #    `graph` - unweighted graph carrying V's
    #
    #    `methods`
    #         list with two functions evaluating different contributions
    #         to the equations of motion:
    #
    #         methods[1](V_1, V_2)
    #                 contribution to the rate due to "coupling"
    #                 between V_1 and V_2
    #
    #         methods[2](V)
    #                 local contribution (anisotropy)
    #            
    #    V(1:|graph|) - current distribution of dynamical variables
    #
    #    Ks - the anisotropy constant
    #
    # OUTPUT:
    #   ΔV(1:|graph|) - array of increments

    out = Ks .* methods[2](V)
    for node in vertices(graph)
        Vnode = V[node]

        for neib in neighbors(graph, node)
            out[node] += methods[1](Vnode, V[neib])
        end
    end
    return out
end

function step_rate(graph::SimpleGraph, method::Function, V::Array, Ks::Float64)
    # Evaluates ΔV for a single step
    # This version presumes that there is a single method for coupling and
    # anisotropy and that there is only easy-axis anisotropy
    #
    # This is a hacky function (see how the anisotropy is evaluated
    # using method(V, -V)). It should be used with care as the
    # implemented idea of a "single method" is flawed. A true single
    # method would imply the presence of a field (a fixed
    # vertex). Here, however, the "single method" is used for
    # anisotropy, which presumes a tensor.
    #
    # INPUT:
    #    graph - unweighted graph carrying V's 
    #    method - method to evaluate different contributions
    #   V(1:|graph|) - current distribution of dynamical variables
    #    Ks - anisotropy constant
    #
    # OUTPUT:
    #   ΔV(1:|graph|) - array of increments

    out = Ks .* method(V, -V)
    for node in vertices(graph)
        Vnode = V[node]

        for neib in neighbors(graph, node)
            out[node] += method(Vnode, V[neib])
        end
    end
    return out
end

function trajectories(graph, methods::Array{Function}, Ks, scale, duration, Vini)
    # Advances the graph ()duration - 1) steps forward taking an array
    # of methods passed further down
    #
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

    tran = 1:(duration - 1)
    for tau in tran
        ΔV = scale .* step_rate(graph, methods, V, Ks)
        V += ΔV
        VFull = [VFull V]
    end

    return VFull
end


function trajectories(graph, method::Function, Ks, scale, duration, Vini)
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

    tran = 1:(duration - 1)
    for tau in tran
        ΔV = scale .* step_rate(graph, method, V, Ks)
        V += ΔV
        VFull = [VFull V]
    end

    return VFull
end

function propagate(graph, method::Function, Ks, scale, duration, Vini)
    # Advances the graph duration - 1 steps forward
    # This is the short version, which returns only the final state vector
    #
    # scale - parameter to tweak the dynamics
    # duration - how many time points to evaluate
    # Vini - the initial conditions
    #
    # OUTPUT:
    #   [V[1] .. V[nv(graph)] at t = duration - 1

    V = Vini

    for tau in 1:(duration - 1)
        ΔV = scale .* step_rate(graph, method, V, Ks)
        V += ΔV
    end

    return V
end

function propagate(model::Model, duration, Vini)
    # Advances the model::Model duration - 1 steps forward
    # Returns only the final state vector
    #
    # model - the model description
    # duration - how many time points to evaluate
    # Vini - the initial conditions
    #
    # OUTPUT:
    #   [V[1] .. V[nv(graph)] at t = duration - 1

    V = Vini
    graph = model.graph
    method = model.method
    Ks = model.Ks
    scale = model.scale

    for tau in 1:(duration - 1)
        ΔV = scale .* step_rate(graph, method, V, Ks)
        V += ΔV
    end

    return V
end

function propagateAdaptively(model::Model, duration, Vini)
    # Advances the model::Model at most duration - 1 steps forward
    # with an adaptive time scale
    #
    #####    NOTE: the current state is under construction!!!     ######
    #
    # Returns only the final state vector
    #
    # model - the model description
    # duration - how many time points to evaluate
    # Vini - the initial conditions
    #
    # OUTPUT:
    #   [V[1] .. V[nv(graph)] at t = duration - 1

    tauconst = 0.6
    cconst = 0.5
    excconst = 10 # how many upscales over the initial one are allowed

    V = Vini
    graph = model.graph
    method = model.method
    scale = model.scale

    curEnergy = energy(graph, method, V)

    exccount = 0 # upscales counter
    tau = 0
    while tau < duration
        exccount = 0 # upscales counter
        grad  = step_rate(graph, method, V, model.Ks)

        ΔV = scale .* grad
        grad2 = sum(grad .* grad)
        bestshift = grad2 * scale * cconst

        candEnergy = energy(graph, method, V + ΔV)

        dEr = abs(candEnergy - curEnergy) / (1 + abs(curEnergy))
        if ( bestshift <= 1e-5)
            debug_msg("Saturation exit at tau = $tau with |dV| = $(sqrt(sum(ΔV .* ΔV))), dEr = $dEr , Enew = $candEnergy, Eold = $curEnergy, grad2 = $grad2, and scale = $scale")
            break
        end

        while candEnergy > curEnergy + bestshift * 0.01 && exccount <= excconst
            scale *= tauconst
            ΔV = scale .* grad
            candEnergy = energy(graph, method, V + ΔV)
            exccount += 1
        end

        # if exccount <= 2
        scale = min(model.scale, 2 * scale)
        # end

        # if candEnergy > curEnergy
        #     exccount += 1
        #     if exccount > excconst
        #         println("Upscales limit is exceeded at tau = $tau")
        #         break
        #     end
        #     scale *= tauconst
        # elseif candEnergy > curEnergy - bestshift
        #     V += ΔV
        # else
        #     while candEnergy < curEnergy - bestshift
        #         exccount -= 1
        #         scale /= tauconst
        #         #println("Upscaling: $scale, $exccount")
        #         ΔV = scale.*grad
        #         candEnergy = energy(graph, method, V + ΔV)
        #         bestshift /= tauconst
        #     end
        #     V += ΔV.*tauconst # we reverse the very last rescaling
        # end

        curEnergy = candEnergy
        tau += 1
    end

    return V
end

function energy(graph, method::Function, V::Array)
    # Evaluates the coupling energy corresponding to V
    # Note: this is essential that this energy evaluates only the coupling energy without
    #       any anisotropic terms
    #
    #####            NOTE: It's broken as of now! Need to pass the correct method
    #
    # INPUT:
    #   method here is the energy of the elementary one-edge graph. Of course, this is
    #          not the same method as in the propagation functions (that woould be the minus
    #          gradient of the method passed to this function)

    en = 0
    for edge in edges(graph)
        en += cosine(V[edge.src] - V[edge.dst])
    end

    return en
end

############################################################
#
### Simulation
#
############################################################

function conf_decay(graph, conf::Array, listlen=3)
    # Evaluates the eigenvalue and eigenvector of the most unstable excitation in configuration conf
    #
    # According to configuration conf, we divide graph into two induced subgraphs (isg) and a bipartite graph (big)
    # so that graph = isg1 + isg2 + big (pluses are set unions)
    # For this separation we evaluate the maximal eigenvalue of L(isg1 + isg2) - L(big), where L is the Laplacian
    #
    # INPUT:
    #   graph
    #   conf - a string with configuration encoded in
    #   listlen - how many eigenvalues to be returned
    #
    # OUTPUT:
    # lambda, v - eigenvalue and eigenvector

    NVert = nv(graph)
    if NVert !== length(conf)
        throw(Error("The configuration legnth does not match the size of the graph"))
    end

    L = laplacian_matrix(graph)

    D = zeros(NVert, 1) # to collect sums of rows

    for (i, j) in zip(findnz(L)...)
        L[i,j] *= conf[i] * conf[j]
        if i != j
            D[i] += L[i,j]
        end
    end

    # Correct diagonal elements of the cut Hessian
    for i in 1:NVert
        L[i,i] = -D[i]
    end

    out = eigs(L, nev=listlen, which=:LR)[1]
    return out[1:listlen]
end

function conf_decay_states(graph, conf::Array, listlen=3)
    # Evaluates the eigenvalue and eigenvector of the most unstable excitation in configuration conf
    #
    # According to configuration conf, we divide graph into two induced subgraphs (isg) and a bipartite graph (big)
    # so that graph = isg1 + isg2 + big (pluses are set unions)
    # For this separation we evaluate the maximal eigenvalue of L(isg1 + isg2) - L(big), where L is the Laplacian
    #
    # INPUT:
    #   graph
    #   conf - a string with configuration encoded in
    #   listlen - how many eigenvalues to be returned
    #
    # OUTPUT:
    # lambda, v - eigenvalue and eigenvector

    NVert = nv(graph)
    if NVert !== length(conf)
        throw(Error("The configuration length does not match the size of the graph"))
    end

    L = laplacian_matrix(graph)

    D = zeros(NVert, 1) # to collect sums of rows
    for (i, j) in zip(findnz(L)...)
        L[i,j] *= conf[i] * conf[j]
        if i != j
            D[i] += L[i,j]
        end
    end

    # Correct diagonal elements of the cut Hessian
    for i in 1:NVert
        L[i,i] = -D[i]
    end

    eigvalue, eigvector = eigs(L, nev=1, which=:LR)
    return (eigvalue, eigvector)
end

function scan_for_best_configuration(model::Model, Vc::Array, domain::Float64, tmax::Integer, Ninitial::Integer)
    # Scans a vicinity of Vc with size given by domain in the "Monte Carlo style"
    #
    # Ninitial number of random initial conditions with individual amplitudes varying between ±domain

    G = model.graph
    Nvert = nv(G)

    (bcut, bconf) = get_best_configuration(G, Vc)

    for i in 1:Ninitial
        local cucu::Integer

        Vi = domain .* get_initial(Nvert, (-1, 1))
        Vi .+= Vc;
        Vi[rand((1:Nvert))] *= -1.1  # a random node is flipped as a perturbation
        # Vi .-= sum(Vi)/Nvert

        # VF = propagateAdaptively(model, tmax, Vi);
        VF = propagate(model, tmax, Vi);

        (cucu, cuco) = get_best_configuration(G, roundup(VF))

        if cut(G, cuco) != cucu
            println("INTERNAL ERROR! Cuts are inconsistent!")
        end

        if cucu > bcut
            # println("Improvement by $(cucu - bcut)")
            bcut = cucu
            bconf = cuco
            Vc = bconf
        end
    end

    return (bcut, bconf)
end

function test_branch(model, Vstart, domain, tmax, Ni, max_depth)
    # Tests branch starting from Vstart
    # INPUT:
    #   model
    #   Vstart
    #   domain      - the size of the vicinity to take the initial states from
    #   tmax        - for how long propagate the particular initial conditions
    #   Ni          - the number of trials
    #   max_depth   - the maximal depth

    local bcut::Integer, bconf::Array{Integer}, stepcount::Integer, nextflag::Bool

    Nvert = nv(model.graph)

    nextflag = true

    (bcut, bconf) = get_best_configuration(model.graph, Vstart)

    Vcentral = Vstart
    stepcount = 0
    while nextflag
        local cucu::Integer

        stepcount += 1

        (cucut, cuconf) = scan_for_best_configuration(model, Vcentral, domain, tmax, Ni)
        message(model, "Prelocal cut = $cucut", 2)
        cuconf = local_search(model.graph, cuconf)
        message(model, "Post-local cut = $(cut(model.graph, cuconf))", 2)

        cuconf = local_twosearch(model.graph, cuconf)
        cucut = cut(model.graph, cuconf)
        message(model, "Post-local-2 cut = $cucut", 2)
        if cucut > bcut
            bcut = cucut
            message(model,
                    "Step $stepcount arrived at $bcut with the displacement $(HammingD(cuconf, bconf))", 2)
            bconf = cuconf
            Vcentral = cuconf

            # domain *= 0.9
            Ni += 10
        else
            nextflag = false
        end

        if stepcount > max_depth
            message(model, "Exceeded maximal depth", 2)
            nextflag = false
        end
    end
    return (bcut, bconf)
end

end # end of module Dice
