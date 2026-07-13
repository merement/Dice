## dyn_anisotropy_model.jl
## 
## Methods for treating anisotropy dynamically
## Currently assumes global dynamical direction
## TODO: 1. Direction and magnitude
##       2. Global and local

# Legacy/experimental definitions moved from Dice.jl
const graph_types = Set([:binary, :weighted])
const isotropy_types = Set([:isotropic, :anisotropic, :dynamic])
const noise_types = Set([:noiseless, :thermal, :dynamic])

struct ModelKind
    graph_type::Symbol
    anisotropy::Symbol
    Noise::Symbol

    function ModelKind()
        new(:binary, :isotropic, :noiseless)
    end

    function ModelKind(graph_t::Symbol)
        graph_t in graph_types ||
            throw(ArgumentError("Invalid graph type: $graph_t"))
        new(graph_t, :isotropic, :noiseless)
    end
end

mutable struct Model
    kind     ::ModelKind
    graph    ::ModelGraph
    coupling ::Function
    scale    ::Float64
    silence  ::Int

    Model() = new(ModelKind())
    Model(x::ModelKind) = new(x)
end

Model(graph::ModelGraph, coupling::Function) =
    begin
        M = graph isa SimpleWeightedGraph ?
            Model(ModelKind(:weighted)) :
            Model(ModelKind(:binary))

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

Model(graph::ModelGraph, coupling::Function, scale::Float64, Ks::Float64) =
    begin
        M = Model(graph, coupling, scale)
        @error "Setting up anisotropy is an error in the present version"

        return M
    end

function step_rate_extended(graph::SimpleGraph, coupling::Function,
                            V::Array,
                            Ks::Float64, anisotropy::Function, R::Float64)
    # Evaluate ΔV and ΔR in the extended system
    #
    # INPUT:
    #    graph - unweighted graph carrying V's 
    #    method - method to evaluate different contributions
    #   V(1:|graph|) - current distribution of dynamical variables
    #    Ks - anisotropy constant
    #
    # OUTPUT:
    #   ΔV(1:|graph|) - array of increments

    #    out = Ks .* method(V, -V) # refactoring to scalars
    out  = zeros(length(V))
    r_out = 0
    for node in vertices(graph)
        Vnode = V[node]
        an = Ks * anisotropy(Vnode - R)
        out[node] = an
        r_out += -an
        for neib in neighbors(graph, node)
            out[node] += coupling(Vnode, V[neib])
        end
    end
    return (out, r_out)
end

function trajectories_extended(graph, duration, scale,
                               pair_method::Function, Vini,
                               Ks, one_method::Function, Rini)
    # Advance the graph `duration - 1` steps forward
    # This is the verbose version, which returns the full dynamics
    #
    #
    # OUTPUT:
    #   (VFull, RFull)
    #           VFull = [V(0) V(1) ... V(duration-1)]
    #           RFull - the final state of the added degree of freedom

    VFull = Vini
    V = Vini
    RFull = [Rini]
    R = Rini

    for _ = 1:(duration - 1)
        (ΔV, ΔR) = step_rate_extended(graph, pair_method, V, Ks,
                                      one_method, R)
         V += scale .* ΔV
        R += scale * ΔR
        VFull = [VFull V]
        RFull = [RFull R]
    end

    return (VFull, RFull)
end

function trajectories_extended(model::Model, duration, Vini, Rini)
    # Advances the graph `(duration - 1)` steps forward
    # This is the verbose version, which returns the full dynamics
    #
    # model - the model description
    # duration - how many time points to evaluate
    # V0 - the initial conditions
    #
    # OUTPUT:
    #   VFull = [V(0) V(1) ... V(duration-1)]

    return trajectories_extended(model.graph, duration, model.scale,
                                 model.coupling,
                        Vini, model.Ks, model.anisotropy, Rini)
end

"""
    propagate_extended(graph, duration, scale, pair_method, Vini,
                              Ks, one_method, Rini)

Propagate extended G + 1 model with the explicit specification.

INPUT:
      graph - the model description
      duration - the number of time points
      scale - the length of the single time step
      pair_method - coupling
      Vini - the initial configuration of G
      Ks - the anisotropy costant
      one_method - the anisotropy method
      Rini - the initial state of the additional degree of freedom

OUTPUT:
      (Vfinal, Rfinal) -
              Vfinal - the final configuration of G
              Rfinal - the final state of the additional degree of freedom
"""
function propagate_extended(graph, duration, scale, pair_method, Vini,
                            Ks, one_method, Rini)

    V = Vini
    R = Rini

    for _ = 1:(duration - 1)
        (ΔV, ΔR) = step_rate_extended(graph, pair_method, V, Ks,
                                      one_method, R)
        V += scale .* ΔV
        R += scale * ΔR
    end
    return (V, R)
end

"""
    propagate_extended(model::Model, duration, Vini, Rini)

Propagate extended G + 1 model.

INPUT:
      model - the model description
      duration - the number of time points
      Vini - the initial configuration of G
      Rini - the initial state of the additional degree of freedom

OUTPUT:
      (Vfinal, Rfinal) -
              Vfinal - the final configuration of G
              Rfinal - the final state of the additional degree of freedom
"""
function propagate_extended(model::Model, duration, Vini, Rini)
    graph = model.graph
    Ks = model.Ks
    pair_method = model.coupling
    one_method = model.anisotropy
    scale = model.scale
    
    return propagate_extended(graph, duration, scale, pair_method, Vini,
                              Ks, one_method, Rini)
end
 
