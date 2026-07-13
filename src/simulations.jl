# simulations.jl
# Simulation functions

function conf_decay(graph, conf::Array, listlen=3)
    # Evaluate the instability of conf
    #
    # INPUT:
    #   graph
    #   conf - a string with configuration encoded in
    #   listlen - how many eigenvalues to be returned
    #
    # OUTPUT:
    # array[1:listlen] of the largest eigenvalues

    NVert = nv(graph)
    if NVert !== length(conf)
        throw("The configuration legnth does not match the size of the graph")
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

# TODO(?) function perturbation_spectrum(graph, conf::Array, listlen=3)
#   listlen - how many eigenvalues to be returned
"""
    conf_decay_states(graph, conf::Array{Int8, 1})

Evaluates the eigenvalue and eigenvector of the most unstable
excitation in configuration conf

INPUT:
    `graph`
    `conf` - a string with configuration encoded in

OUTPUT:
    lambda, v - eigenvalue and eigenvector

"""
function conf_decay_states(graph, conf::Array{Int8, 1})

    NVert = nv(graph)
    if NVert !== length(conf)
        throw("The configuration length does not match the size of the graph")
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

function scan_for_best_configuration(model::Model, Vc::Array,
                                     domain::Float64, tmax::Integer,
                                     Ninitial::Integer)
# Scan a vicinity of Vc with size given by domain in the "Monte Carlo style"
#
# Ninitial number of random initial conditions with individual amplitudes
# varying between ±domain

    G = model.graph
    (bcut, bconf) = get_best_configuration(G, Vc)
    Nvert = nv(G)
    for _ = 1:Ninitial
#        local cucu::Integer
        Vi = domain .* get_random_cube(Nvert, (-1, 1))
        Vi .+= Vc;
        Vi[rand((1:Nvert))] *= -1.1  # flip a random node as a perturbation
        # Vi .-= sum(Vi)/Nvert

        # VF = propagateAdaptively(model, tmax, Vi);
        VF = propagate(model, tmax, Vi);

        (cucu, cuco) = get_best_configuration(G, roundup(VF))

        if cucu > bcut
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

#    local bcut::Integer, bconf::Array{Integer}, stepcount::Integer, nextflag::Bool

    (bcut, bconf) = get_best_configuration(model.graph, Vstart)

    Vcentral = Vstart
    stepcount = 0
    nextflag = true
    while nextflag
#        local cucu::Integer
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
                    "Step $stepcount arrived at $bcut with the displacement $(hamming_distance(cuconf, bconf))", 2)
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

