# Script deprecated.jl

# Collection of deprecated methods in a random order
# Can be moved to obsoletes.jl unpredictably, do not use!

function get_initial(Nvert::Int, (vmin, vmax))
    # Generate random vector with Nvert components uniformly
    #  distributed in the (vmin, vmax) interval
    # NOTE: OBSOLETE, replaced by get_random_sphere or get_random_cube
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

function realign_2(conf::Hybrid, r::Float64 = 0.0)
    # Changes the reference point for the separated representation by `r`
    # according to xi - r = sigma(r) + X(r)
    # INPUT & OUTPUT:
    #     conf = (sigma, X)
    #
    # NOTE: Obsolete, replaced by realign_hybrid
    V = Dice.hybrid_to_cont(conf)
    return Dice.cont_to_hybrid(V, r)
end

function update_2!(spins::SpinConf, xs::FVector, dx::FVector)
    # Tracing version
    # Update the continuous component (xs) by dx using the wrapping rule
    # The spin part is assumed to be Int8-array
    # Return the number of flips (tracing)
    # NOTE: deprecated in favor of "hybrid" methods
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
    for _ = 1:(tmax-1)
        #        DX = step_rate_hybrid(graph, mtd, S, X, Ns, noise).*scale
        DX = step_rate_hybrid(graph, mtd, S, X) .* scale
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
    for _ = 1:(tmax-1)
        #        DX = step_rate_hybrid(graph, mtd, S, X, Ns, noise).*scale
        DX = step_rate_hybrid(graph, mtd, S, X) .* scale
        update_2!(S, X, DX)
        XFull = [XFull X]
        SFull = [SFull S]
    end
    return (SFull, XFull)
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
function separate(Vinp::Array{Float64,1}, r=0.0)
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
function combine(s::Array{Int8}, x::Array{Float64}, r=0.0)
    return x .+ s .+ r
end

# Backward compatibility wrappers for renamed functions
export
    saveMTXGraph,
    loadMTXGraph,
    loadMTXAdjacency,
    HammingD,
    EuclidD,
    noiseUniform,
    noiseNormal

function saveMTXGraph(G::ModelGraph, filename::AbstractString)
    @warn "Function `saveMTXGraph` is deprecated, use `save_mtx_graph` instead"
    save_mtx_graph(G, filename)
end

function loadMTXGraph(filename::AbstractString, min_weight::Float64=weight_eps)::ModelGraph
    @warn "Function `loadMTXGraph` is deprecated, use `load_mtx_graph` instead"
    load_mtx_graph(filename, min_weight)
end

function loadMTXAdjacency(filename::AbstractString, min_weight::Float64=weight_eps)::Matrix
    @warn "Function `loadMTXAdjacency` is deprecated, use `load_mtx_adjacency` instead"
    load_mtx_adjacency(filename, min_weight)
end

function HammingD(s1::AbstractVector, s2::AbstractVector)::Int
    @warn "Function `HammingD` is deprecated, use `hamming_distance` instead"
    hamming_distance(s1, s2)
end

function EuclidD(V1::AbstractVector{<:Real}, V2::AbstractVector{<:Real})::Float64
    @warn "Function `EuclidD` is deprecated, use `euclid_distance` instead"
    euclid_distance(V1, V2)
end

function noiseUniform(L::Int)::Array{Float64, 1}
    @warn "Function `noiseUniform` is deprecated, use `noise_uniform` instead"
    noise_uniform(L)
end

# Candidates to deprecate (Moved from src/Dice.jl)
export
    get_random_cut,
    get_rate,
    cut_2

function get_random_cut(graph, V, trials=10)
    @warn "Function `get_random_cut` is deprecated and will be removed in a future version"
    (becu, _, _) = get_random_rounding(graph, V, trials)
    return becu
end

function get_rate(VFull)
    @warn "Function `get_rate` is deprecated and will be removed in a future version"
    out = [sum((VFull[:, i+1] - VFull[:, i]) .^ 2) for i in 1:(size(VFull)[2]-1)]
    return sqrt.([out; out[end]])
end

function cut_2(model, s::SpinConf, x::FVector)
    @warn "Function `cut_2` is deprecated and will be removed in a future version"
    phix = 0
    graph = model.graph
    for edge in edges(graph)
        phix += s[edge.src] * s[edge.dst] * abs(x[edge.src] - x[edge.dst]) / 2
    end
    return Dice.cut(graph, s) + phix
end

function cut_2(model, distr::Hybrid)
    return cut_2(model, distr[1], distr[2])
end

