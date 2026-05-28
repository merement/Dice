# project: logic
# module RCAdder.jl

# General functions for the relaxed-spin implementation of RCA
# NOTE: The current version presumes a specific relation between the
# logical bits and the spin representation as defined in RSLogic
#

module RCAdder

# see "more general" in RSLogic
# node convention: (z, x, y, s, c)
# weights (zx, zy, zs, zc, xy, xs, xc, ys, yc, sc)
# full adder weights (from [1]) (used from RSLogic)
# const FAA_weights_0 = [1, 1, -1, -2, 1, -1, -2 * (1 + 0.0), -1, -2 * (1 + 0.0), 2]
# full adder weights (updated)
const FAA_weights = [1, 1, -1, -2, 1, -1, -2 * (1 + 0.1), -1, -2 * (1 + 0.1), 2]

using RSLogic

using Dice, Graphs, SimpleWeightedGraphs
using RSGraphics
using Primes

# constants
const MAX_X = 1.0 # Hard constant - the right end of the continuous component
const MIN_X = -1.0 # Hard constant - the left end of the continuous component
const BASE_X = 0.0 # Soft constant - the x-coordinate of "fixed" elements

# The displaced spins sector
const SECTOR_B = -0.5 # Position of the (x, y) pair (the bottom of the sector)
const SECTOR_T = 0.5 # The upper boundary of the sector

const DELTA_V = 0.1 # half-size of strong cluster (tolerance)

const IVector = Vector{Int}
const I8Vector = Vector{Int8}
const FVector = Vector{Float64}

# Spin and decimal representations of the result
# (input_1, input_2, sum, carry-in, carry-out)
const Result = Tuple{I8Vector, I8Vector, I8Vector, Int8, Int8}
const ResultDec = Tuple{Int, Int, Int, Int, Int}


export
    generate_rcadder_graph,
    addition_to_spin_state

function prime_list(num ::Integer) ::IVector
    list_p ::IVector = []
    value = 2
    for _ in 1:num
        push!(list_p, value)
        value = Primes.nextprime(value, 2)
    end
    return list_p
end

function get_independent_weights(Nbit ::Int)
    p_list = prime_list(Nbit)
    l_list = sqrt.(p_list)
    return l_list ./= round.(l_list)
end

function get_rcadder_graph(Nbit ::Int, FA_weights ::FVector, l_list ::FVector)
    num_vertices = 5 # number of vertices in FA
    adder_vertices = 4 * Nbit + 1
    add_graph = SimpleWeightedGraph(adder_vertices)
    start = 1
    for fa in 1:Nbit
        w_position = 1
        for node in start:(start + num_vertices - 2)
            for recipient in (node + 1):(start + num_vertices - 1)
                add_edge!(add_graph,
                          node, recipient,
                          l_list[fa] * FA_weights[w_position])
                w_position += 1
            end
        end
        start += 4
    end
    return add_graph
end


# NOTE: it uses weights hardcoded here in the module
function get_rcadder_weighted_graph(Nbit ::Int, l_list ::FVector)
    num_vertices = 5 # number of vertices in FA
    adder_vertices = 4 * Nbit + 1
    add_graph = SimpleWeightedGraph(adder_vertices)
    start = 1
    for fa in 1:Nbit
        w_position = 1
        for node in start:(start + num_vertices - 2)
            for recipient in (node + 1):(start + num_vertices - 1)
                add_edge!(add_graph,
                          node, recipient,
                          l_list[fa] * FAA_weights[w_position])
                w_position += 1
            end
        end
        start += 4
    end
    return add_graph
end

function get_rcadder_weighted_graph(Nbit ::Int)
    return get_rcadder_weighted_graph(Nbit, get_independent_weights(Nbit))
end

function generate_rcadder_graph(Nbit ::Int)
    return get_rcadder_weighted_graph(Nbit, ones(Nbit))
end

# Generate Nbit adder graph with all "default" weights
function generate_rcadder_0_graph(Nbit ::Int)
    return get_rcadder_graph(Nbit,
                             1.0 .* RSLogic.generate_gate(:FA)[2],
                             ones(Nbit))
end

# Convert a _spin_ form of an addition problem (`add1 + add2 + cin`) with
# the optional carry-in `cin` to the spin representation for the chain of
# hard connected 5-spin FAs:
# (cin, add1, add2, sum, carry, add1, add2, sum, ...).
# Free spins are set to 0. Their initialization is delegated to further
# processing.
function addition_to_spin_state(add1 ::I8Vector,
                                add2 ::I8Vector,
                                cin = -1) ::I8Vector
    num_bits = length(add1)
    add2 = add2[1:num_bits]
    spin_state ::I8Vector = []
    
    push!(spin_state, cin)
    for k in 1:num_bits
        push!(spin_state, add1[k])
        push!(spin_state, add2[k])
        push!(spin_state, Int8(0)) # sum
        push!(spin_state, Int8(0)) # carry
    end        
    return spin_state     
end

# convert an _int_ form of an addition problem
function addition_to_spin_state(add1 ::Int, add2 ::Int, cin = 0) ::I8Vector
    a, b = ints_to_spin_arrays(add1, add2)
    cin = cin == 0 ? -1 : 1
    return addition_to_spin_state(a, b, cin)
end

# Take a spin representation of the addition problem (for instance, produced
# by `addition_to_spin_state`) and produce the randomized initial state of
# the network
#
# NOTE: this function spreads displaced bits inside a preallocated sector.
# For distributing over a list of specific positions, use
# get_initial_rca_rnd_marks!
#
# `inv_seq` specifies the sequence of displaced r-spins, 0-based enumerated
# `0` corresponds to the carry-in
# 'k' corresponds to the k-th bit pair
function get_initial_rca_rnd!(add_state ::I8Vector,
                              inv_seq ::IVector) ::Tuple{Tuple{I8Vector,
                                                               FVector},
                                                         IVector}
    num_spins = length(add_state)
    num_preblock = num_spins - 1
    mod(num_preblock, 4) != 0 && error("Invalid spin state is provided")
    # 0-based number of input bits/bit-pairs
    num_bits = div(num_preblock, 4)
    
    # First, the discrete component, then, we'll deal with variations
    free_spins = findall(x -> x == 0, add_state)
    num_free = length(free_spins)
    nonfree = setdiff(1:num_spins, free_spins)
    # randomize free discrete components
    add_state[free_spins] = rand(Int8[-1, 1], num_free)
    
    x_comp = SECTOR_B .* ones(num_spins) # default to SECTOR_B
    
    num_invs = length(inv_seq)
    if num_invs == 0
        x_comp[free_spins] = 2 .* rand(num_free) .- 1 # frees are randomized
        return ((add_state, x_comp), nonfree)
    end
    x_comp[free_spins] = Dice.get_random_cube(length(free_spins), 1.99)

    (minimum(inv_seq) < 0 || maximum(inv_seq) > num_bits) &&
        error("Invalid sequence of inversions: Out of bounds")    
    num_invs != length(unique(inv_seq)) &&
        error("Invalid sequence of inversions: Repeated elements")

    delta_x = (SECTOR_T - SECTOR_B) / num_invs
    
    for (inv_ind, inv_bits) in enumerate(inv_seq)
        spin_pos = inv_bits == 0 ?
            [1] :
            bit_pos_to_spin_pos(inv_bits)
        complex_x = SECTOR_B + inv_ind * delta_x
        x_comp[spin_pos] .= complex_x
    end
    state = (add_state, x_comp)
    return (state, nonfree)
end

function get_initial_rca_rnd(add_state ::I8Vector,
                             inv_seq ::IVector) ::Tuple{Tuple{I8Vector,
                                                              FVector},
                                                        IVector}
    problem = copy(add_state)
    return get_initial_rca_rnd!(problem, inv_seq)
end

function get_initial_rca_rnd(add1,
                             add2,
                             inv_seq ::IVector) ::Tuple{Tuple{I8Vector,
                                                               FVector},
                                                         IVector}
    return get_initial_rca_rnd!(addition_to_spin_state(add1, add2), inv_seq)
end

# marks - the list of coordinates to place (base and displaced) fixed spins
# The list must contain at least one (base) element
# inv_seq - the list of invertable bits pairs (or carry-in). Convention: the
# list has elements in the order of their inversions while rotating the
# phase circle counterclock wise
# Example: inv_list = [a, b, c] then X_a > X_b > X_c
# with marks = [X_base, ..., X_c, X_b, X_a]
function get_initial_rca_rnd_marks!(add_state ::I8Vector,
                                   inv_seq ::IVector,
                                   marks ::FVector) ::Tuple{Tuple{I8Vector,
                                                                  FVector},
                                                            IVector,
                                                            FVector}
    length(marks) < 1 &&
        error("Invalid marks list: must contain at least the base element")
    
    num_spins = length(add_state)
    num_preblock = num_spins - 1
    mod(num_preblock, 4) != 0 && error("Invalid spin state")
    # 0-based number of input bits/bit-pairs
    num_bits = div(num_preblock, 4)

    num_invs = length(inv_seq)
    if num_invs + 1 != length(marks)
        @show (num_invs, marks)
        error("Logical and physical positions of displaced spins don't match")
    end
    
    # First, the discrete component
    free_spins = findall(x -> x == 0, add_state)
    num_free = length(free_spins)
    nonfree = setdiff(1:num_spins, free_spins)
    add_state[free_spins] = rand(Int8[-1, 1], num_free)

    # X-component: free randomized, non-free set to marks
    x_comp = marks[1] .* ones(num_spins) # default - nondisplaced
    x_comp[free_spins] = get_random_x(num_free) 
    if num_invs == 0
        return ((add_state, x_comp), nonfree, marks)
    end
    
    (minimum(inv_seq) < 0 || maximum(inv_seq) > num_bits) &&
        error("Invalid sequence of inversions: Out of bounds")    
    num_invs != length(unique(inv_seq)) &&
        error("Invalid sequence of inversions: Repeated elements")

    for (inv_ind, inv_bits) in enumerate(inv_seq)
        spin_pos = inv_bits == 0 ?
            [1] :
            bit_pos_to_spin_pos(inv_bits)
        x_comp[spin_pos] .= marks[end - inv_ind + 1]
    end
    state = (add_state, x_comp)
    return (state, nonfree, marks)
end

function get_initial_rca_rnd_marks!(add_state ::I8Vector,
                                    inv_seq ::IVector) ::Tuple{Tuple{I8Vector,
                                                                     FVector},
                                                               IVector,
                                                               FVector}
    return get_initial_rca_rnd_marks!(add_state,
                                     inv_seq,
                                     get_marks(inv_seq))
end

function get_initial_rca_rnd_marks(add_state ::I8Vector,
                                   inv_seq ::IVector) ::Tuple{Tuple{I8Vector,
                                                                    FVector},
                                                              IVector,
                                                              FVector}
    problem = copy(add_state)
    return get_initial_rca_rnd_marks!(problem, inv_seq)
end

function get_initial_rca_rnd_marks(add1,
                                   add2,
                                   inv_seq ::IVector) ::Tuple{Tuple{I8Vector,
                                                                    FVector},
                                                              IVector,
                                                              FVector}
    return get_initial_rca_rnd_marks!(addition_to_spin_state(add1, add2), inv_seq)
end


# extract from a (binary) spin state operands, the sum, the carry-in
# and the carry-out
function spin_state_to_result(state ::I8Vector) ::Result
    num_spins = length(state)
    num_preblock = num_spins - 1
    mod(num_preblock, 4) != 0 && error("Invalid spin state is provided")
    num_bits = div(num_preblock, 4)
    
    cin = state[1]
    a_arr = state[2:4:end]
    b_arr = state[3:4:end]
    s_arr = state[4:4:end]
    cout = state[end]
    return (a_arr, b_arr, s_arr, cin, cout)
end

# extract from a (binary) spin state operands, the sum, and the carry-out
# return the results in the decimal representation
function spin_state_to_result_dec(state ::I8Vector) ::Result
    a, b, s, cin, cout = spin_state_to_result(state)
    return (spins_to_int.((a, b, s)),
            cin,
            cout)
end

function state_to_result(rs_state ::Tuple{I8Vector, FVector},
                         inv_seq ::IVector) ::Vector{Result}
    out_res ::Vector{Result} = []
    
    spin_state = rs_state[1]
    push!(out_res, spin_state_to_result(spin_state))

    num_invs = length(inv_seq)
    if num_invs == 0
        # we always have at least global inversion
        push!(out_res, spin_state_to_result(-spin_state))
        return out_res
    end
    
    num_spins = length(spin_state)
    num_bits = div(num_spins - 1, 4)
    (minimum(inv_seq) < 0 || maximum(inv_seq) > num_bits) &&
        error("Invalid sequence of inversions: Out of bounds")    
    num_invs != length(unique(inv_seq)) &&
        error("Invalid sequence of inversions: Repeated elements")
    num_invs > num_spins &&
        error("Invalid sequence of inversions: Too many inversions")
    
    delta_x = (SECTOR_T - SECTOR_B) / num_invs
    
    # collection of nodes associated with displaced bits
    members ::Vector{IVector} = []

    # the list of bounds identifying inverstion intervals
    bounds = get_sectors_from_marks(inv_seq)
    for i in 1:(num_invs + 1)
        push!(members, findall(x -> bounds[i + 1] <= x < bounds[i], rs_state[2]))
    end
    
    # check that all nodes are accounted for
    rem_nodes = 1:num_spins
    for complex in members
        rem_nodes = setdiff(rem_nodes, complex)
    end
    if length(rem_nodes) > 0 
        @show bounds, members, rem_nodes
        error("Not all spins got covered in inversions")
    end
    
    # sequential inversions
    for complex in members
        spin_state[complex] .*= -1
        push!(out_res, spin_state_to_result(spin_state))
    end
    return out_res        
end

function state_to_result_marks(rs_state ::Tuple{I8Vector, FVector},
                               inv_seq ::IVector,
                               marks ::FVector) ::Vector{Result}
    # NOTE: there's apparently a universal validation block
    num_invs = length(inv_seq)
    num_invs + 1 != length(marks) &&
        error("Logical and physical positions of displaced spins don't match")
    spin_state = rs_state[1]
    num_spins = length(spin_state)
    num_bits = div(num_spins - 1, 4)
    
    out_res ::Vector{Result} = []
    push!(out_res, spin_state_to_result(spin_state))

    if num_invs == 0
        # we always have at least global inversion
        push!(out_res, spin_state_to_result(-spin_state))
        return out_res
    end    

    (minimum(inv_seq) < 0 || maximum(inv_seq) > num_bits) &&
        error("Invalid sequence of inversions: Out of bounds")    
    num_invs != length(unique(inv_seq)) &&
        error("Invalid sequence of inversions: Repeated elements")
    num_invs > num_spins &&
        error("Invalid sequence of inversions: Too many inversions")

    # collection of nodes associated with displaced bits
    members ::Vector{IVector} = []

    # the list of bounds identifying inverstion intervals
    bounds = get_sectors_from_marks(marks)
    for i in 1:(num_invs + 1)
        complex = findall(x -> bounds[i + 1] <= x < bounds[i], rs_state[2])
        push!(members, complex)
    end
    
    rem_nodes = 1:num_spins
    for complex in members
        rem_nodes = setdiff(rem_nodes, complex)
    end

    # If anything's left, we've got a problem
    if length(rem_nodes) > 0 
        @show rs_state
        @show bounds
        @show members, rem_nodes
        @show rs_state[2][rem_nodes]
        error("Not all spins got covered in inversions")
    end
    
    # sequential inversions
    for complex in members
        spin_state[complex] .*= -1
        push!(out_res, spin_state_to_result(spin_state))
    end
    return out_res        
end

# extract results from the spin form of the sum state
# This version does not consider the last spin as the carry-out
# NOTE: is it really needed? (Incomplete)
function spin_sum_state_to_result_marks(rs_state ::Tuple{I8Vector, FVector},
                                    inv_seq ::IVector,
                                    marks ::FVector) ::Vector{Result}
    # NOTE: there's apparently a universal validation block
    num_invs = length(inv_seq)
    num_invs + 1 != length(marks) &&
        error("Logical and physical positions of displaced spins don't match")
    spin_state = rs_state[1]
    num_spins = length(spin_state)
    num_bits = div(num_spins - 1, 4)
    
    out_res ::Vector{Result} = []
    push!(out_res, spin_state_to_result(spin_state))

    if num_invs == 0
        # we always have at least global inversion
        push!(out_res, spin_state_to_result(-spin_state))
        return out_res
    end    

    (minimum(inv_seq) < 0 || maximum(inv_seq) > num_bits) &&
        error("Invalid sequence of inversions: Out of bounds")    
    num_invs != length(unique(inv_seq)) &&
        error("Invalid sequence of inversions: Repeated elements")
    num_invs > num_spins &&
        error("Invalid sequence of inversions: Too many inversions")

    # collection of nodes associated with displaced bits
    members ::Vector{IVector} = []

    # the list of bounds identifying inverstion intervals
    bounds = get_sectors_from_marks(marks)
    for i in 1:(num_invs + 1)
        complex = findall(x -> bounds[i + 1] <= x < bounds[i], rs_state[2])
        push!(members, complex)
    end
    
    rem_nodes = 1:num_spins
    for complex in members
        rem_nodes = setdiff(rem_nodes, complex)
    end

    # If anything's left, we've got a problem
    if length(rem_nodes) > 0 
        @show rs_state
        @show bounds
        @show members, rem_nodes
        @show rs_state[2][rem_nodes]
        error("Not all spins got covered in inversions")
    end
    
    # sequential inversions
    for complex in members
        spin_state[complex] .*= -1
        push!(out_res, spin_state_to_result(spin_state))
    end
    return out_res        
end


function result_to_dec(res ::Result)  ::ResultDec
    return ( (x -> bits_to_int(Int.(reverse(div.(x .+ 1, 2))))).(res[1:3])...,
             res[4],
             res[5])
end

function state_to_result_dec(rs_state ::Tuple{I8Vector, FVector},
                             inv_seq ::IVector) ::Vector{ResultDec}
    result = state_to_result(rs_state, inv_seq)
    return result_to_dec.(result)
end

function state_to_result_marks_dec(rs_state ::Tuple{I8Vector, FVector},
                                  inv_seq ::IVector,
                                  marks ::FVector) ::Vector{ResultDec}
    result = state_to_result_marks(rs_state, inv_seq, marks)
    return result_to_dec.(result)
end

function get_sum_state(rs_state ::Dice.Hybrid) ::Dice.Hybrid
    outcome_s = rs_state[1][4:4:end]
    outcome_x = rs_state[2][4:4:end]
    return (outcome_s, outcome_x)
end

#########
##
## Reporting functions
##
#########

# return correct/incorrect, obtained result, expected result
function check_result(res ::Result) ::Tuple{Bool, ResultDec, ResultDec}
    num_bits = length(res[1])
    s_max = 2^num_bits - 1

    input = spins_to_int.(res[1:2])
    sum_res = spins_to_int(res[3])
    cin_inp = res[4]
    cout_res = res[5]
    rep_res = (input[1], input[2], sum_res, cin_inp, cout_res)

    cin_val = cin_inp == -1 ? 0 : 1
    sum_exp = sum(input) + cin_val
    cout_exp = -1
    if sum_exp > s_max
        sum_exp -= (s_max + 1)
        cout_exp = 1
    end
    rep_exp = (input[1], input[2], sum_exp, cin_inp, cout_exp)
    status = (sum_res == sum_exp) && (cout_res == cout_exp)
    return status, rep_res, rep_exp
end

function is_result_correct(res ::Result) ::Bool
    return check_result(res)[1]
end

function report_result(res ::Result) ::Bool
    check = check_result(res)
    num_bits = length(res[1])
    over = 2^num_bits

    out_str = "Tested $(check[3][1]) + $(check[3][2]) + <$(check[3][4])> "
    println(out_str)
    out_str = "\tExpected $(check[3][3]) + <$(check[3][5])> "
    if check[3][5] > 0
        add_str = "($(over + check[3][3]))"
        out_str *= add_str
    end
    println(out_str)
    out_str = "\tObtained $(check[2][3]) + <$(check[2][5])> "
    if check[2][5] > 0
        add_str = "($(over + check[2][3]))"
        out_str *= add_str
    end
    println(out_str)
    if check[1]
        out_str = "Test passed"
    else
        out_str = "Test failed"
    end
    println(out_str)
    
    return check[1]
end

#########
##
## Service functions
##
#########

# Generate a uniform distribution U{[a,b] ∪ [c, d]}
# over two non-overlapping intervals
# [a, b] < [c, d]
function rand_two_intervals(a, b, c, d, n)
    !(a < b && b < c && c < d) &&
        error("Invalid intervals: the boundaries must satisfy a < b < c < d")
    len1 = b - a
    len2 = d - c
    total_len = len1 + len2

    out = rand(n) .* total_len
    part1 = findall(x -> x < len1, out)
    part2 = setdiff(1:n, part1)
    out[part1] .+= a
    out[part2] .+= c - len1
    return out
end

# Convert 1-based position of a bit pair in the pair of operands to a pair of
# numbers of spins in the spin representation.
# `k=0` is a special case corresponding to carry-in
function bit_pos_to_spin_pos(k ::Int) ::IVector
    k == 0 && return [1]
    
    k -= 1 # now, it's 0-based
    x_pos = 1 + 4 * k + 1
    y_pos = x_pos + 1
    return [x_pos, y_pos]
end

function int_to_bit_array(a::Int)
    bits_a = string(a, base=2)
    return [parse(Int, c) for c in bits_a]
end

function int_to_spin_array(a ::Int) ::I8Vector
    arr_a = int_to_bit_array(a)
    arr_a = 2 .* Int8.(arr_a) .- 1
    return reverse(arr_a)
end

function ints_to_bit_arrays(a::Int, b::Int)
    # Get binary representations
    bits_a = string(a, base=2)
    bits_b = string(b, base=2)
    
    L = max(length(bits_a), length(bits_b))    
    # Pad with leading zeros
    padded_a = lpad(bits_a, L, '0')
    padded_b = lpad(bits_b, L, '0')
    
    arr_a = [parse(Int, c) for c in padded_a]
    arr_b = [parse(Int, c) for c in padded_b]
    
    return arr_a, arr_b
end

function ints_to_spin_arrays(a ::Int, b ::Int) ::Tuple{I8Vector, I8Vector}
    arr_a, arr_b = ints_to_bit_arrays(a, b)
    arr_a = 2 .* Int8.(arr_a) .- 1
    arr_b = 2 .* Int8.(arr_b) .- 1
    return reverse(arr_a), reverse(arr_b)
end

function bits_to_int(bits::IVector)
    parse(Int, join(bits); base=2)
end

function bits_to_int(bits::I8Vector)
    parse(Int, join(bits); base=2)
end

function spins_to_int(spins::I8Vector)
    bits = reverse(div.(spins .+ 1, 2))
    parse(Int, join(bits); base=2)
end

"""
    invert_bit(n::Integer, k::Integer) -> Integer

Return `n` with bit `k` inverted (toggled), where `k = 1` is
the least-significant bit.
(So `k = 2` toggles the next bit, etc.)

NOTE: This function works with 64-bit signed integers.
`k > 64` will be silently ignored.
"""
function invert_bit(n::Integer, k::Integer)
    if k < 1
        throw(ArgumentError("bit index k must be >= 1 (LSB is k=1)"))
    end

    mask = one(n) << (k - 1)
    return Base.xor(n, mask)
end

"""
    invert_bits(n::Integer, ks::Vector{Integer}) -> Vector{Integer}

Given an integer `n` and a list of bit positions `ks`, return a vector
where the i-th element is `n` with bit `ks[i]` inverted. Bit positions
are 1-based (LSB is 1).

NOTE: This function works with 64-bit signed integers.
`k > 64` will be silently ignored.
"""
function invert_bits(n::Integer, ks::Vector{Integer})
    out = Vector{Integer}(undef, length(ks))
    for (i, k) in enumerate(ks)
        if k < 1
            throw(ArgumentError("bit index k must be >= 1 (LSB is k=1)"))
        end
        mask = one(n) << (k - 1)
        out[i] = Base.xor(n, mask)
    end
    return out
end

"""
    build_chain(t::Tuple{Integer,Integer,Integer}, ks::Vector{Integer}) -> Vector{Tuple{Integer,Integer,Integer}}

`t` is (c, x, y) with c in {0,1}. Bit positions are 1-based (LSB is 1).
Output[1] == t.
For l = 1:length(ks):
- if ks[l] == 0: toggle c
- if ks[l] > 0: toggle bit ks[l] in both x and y

Validity checks:
- c must be 0 or 1
- all entries of ks must be >= 0
- all entries of ks must be pairwise distinct

NOTE: This function works with 64-bit signed integers.
`k > 64` will be silently ignored.
"""
function build_chain(t::Tuple{Int,Int,Int},
                     ks::Vector{Int}) ::Vector{Tuple{Int,Int,Int}}
    c, x, y = t

    if !(c == 0 || c == 1)
        throw(ArgumentError("c must be 0 or 1"))
    end

    if !isempty(ks) && minimum(ks) < 0
        throw(ArgumentError("ks entries must be >= 0"))
    end

    if length(unique(ks)) != length(ks)
        throw(ArgumentError("all elements in ks must be different"))
    end

    out = Vector{Tuple{Integer,Integer,Integer}}(undef, length(ks) + 1)
    out[1] = (c, x, y)

    for (l, k) in enumerate(ks)
        c_prev, x_prev, y_prev = out[l]

        if k == 0
            c_new = Base.xor(c_prev, 1)
            out[l + 1] = (c_new, x_prev, y_prev)
        else
            mask = one(x_prev) << (k - 1)
            x_new = Base.xor(x_prev, mask)
            y_new = Base.xor(y_prev, mask)
            out[l + 1] = (c_prev, x_new, y_new)
        end
    end

    return out
end

function get_random_x(num_spins ::Int)
    return 2 .* rand(num_spins) .- 1
end

# generate the list of the inversion intervals boundaries 
function get_sectors_from_marks(marks ::FVector)
    middles = [(marks[i] + marks[i + 1])/2 for i in 1:(length(marks) - 1)]
    return reverse(vcat([-1.0], middles, [1.0]))
end

function get_sectors_from_marks(inv_seq ::IVector)
    return get_sectors_from_marks(get_marks(inv_seq))
end

# return positions of fixed spins when `num_sec` spins are displaced
# 2025-12-21: checked and confirmed
function get_marks(num_sec ::Int) ::FVector
    num_sec <= -1 &&
        error("Invalid number of displaced spins (must be positive)")

    num_points = num_sec + 1
    dx_arc = 2 / num_points
    return collect(((1:num_points) .- 0.5) .* dx_arc .- 1)
end

function get_marks(inv_seq ::IVector) ::FVector
    return get_marks(length(inv_seq))
end

function get_rnd_number(nbits ::Int, flag = false) ::Vector{Int}
    out ::Vector{Int} = rand(Bool, nbits)
    if flag
        out[1] = 1
    end
    return out
end

end # module ends here

# References:
#
# [1] S. Tsukiyama et al., “Designing unit Ising models for logic gate
# simulation through integer linear programming,” International Journal of
# Parallel, Emergent and Distributed Systems, vol. 40, no. 2, pp. 116–136,
# Mar. 2025, doi: 10.1080/17445760.2024.2438043.
