# project: logic
# module RSLogic.jl

# Coommon definitions for relaxed-spin implmenentations of logical gates

module RSLogic

export
    generate_gate

# GATES ::Dict(Symbol, Tuple{Int, Vector{Int}}) =
#     {:AND => (4, [-1, -1, 2, 1, -2, -2])}

# GATES[:FA] = (5, [1, 1, -1, -2, 1, -1, -2, -1, -2, 2])

function generate_gate(type ::Symbol) ::Tuple{Int, Vector{Int}}
    if type == :AND
        # 4-spin (including the aux-spin) implementation of AND gate
        # x ∧ y = a
        # Source: gate H_{AND} from Table 2 in [1]
        # node convention: (f, x, y, a) (f - aux-spin)
        # weights are the upper triangle of the adjacency matrix
        # weights: (xf, yf, af, xy, xa, ya)
        num_vert = 4
        g_weights = [-1, -1, 2, 1, -2, -2]
        return (num_vert, g_weights)
    end

    if type == :FA
        # 5-spin implementation of full adder gate
        # z + x + y = s + 2 * c
        # z - carry-in
        # (x, y) - input bits
        # s - sum
        # c - carry-out
        # Source: gate H_{FA}(x, y, z, s, c) from Table 2 in [1]
        # node convention: (z, x, y, s, c)
        # weights (zx, zy, zs, zc, xy, xs, xc, ys, yc, sc)
        FAA_weights = [1, 1, -1, -2, 1, -1, -2, -1, -2, 2]
        num_vert = 5
        return (num_vert, FAA_weights)
    end

    if type == :FAFixed
        # WARN: doesn't work (internal inconsistency between the overcharged
        # and neutral cases)
        # 5-spin implementation of full adder gate
        # z + x + y = s + 2 * c
        # z - carry-in
        # (x, y) - input bits
        # s - sum
        # c - carry-out
        # Source: gate H_{FA}(x, y, z, s, c) from Table 2 in [1]
        # node convention: (z, x, y, s, c)
        # weights (zx, zy, zs, zc, xy, xs, xc, ys, yc, sc)
        FAA_weights = [0, 0, -3, -3, 0, -1/2, -4, -1/2, -4, 1]
        num_vert = 5
        return (num_vert, FAA_weights)
    end

    if type == :FAaux
        # 6-spin implementation of full adder gate
        # z + x + y = s + 2 * c
        # z - carry-in
        # (x, y) - input bits
        # s - sum
        # c - carry-out
        # u - free aux-spin
        # Source: gate H_{FA}^{zu}(x, y, z, s, c, u) from Table 2 in [1]
        # HFA = c s - s u - c x - u x - c y - s y + u y - c z - u z + x z
        # node convention: (z, x, y, u, s, c)
        # weights (zx, zy, zu, zs, zc, xy, xu, xs, xc, yu, ys, yc, us, uc, sc)
        FAA_weights = [1, 0, -1, 0, -1, 0, -1, 0, -1, 1, -1, -1, -1, 0, 1]
        num_vert = 6
        return (num_vert, FAA_weights)
    end

    error("Unsupported gate $type")
end

# function generate_gate(filename ::String)
# end

end

# References:
#
# [1] S. Tsukiyama et al., “Designing unit Ising models for logic gate
# simulation through integer linear programming,” International Journal of
# Parallel, Emergent and Distributed Systems, vol. 40, no. 2, pp. 116–136,
# Mar. 2025, doi: 10.1080/17445760.2024.2438043.
