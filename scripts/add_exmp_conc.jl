# Concurrent circuit
# Full adders are processed simultaneously
#
# The dependence of the probability of obtaining correct result
# for each displacement group on the duration of each stage

include("intro.jl")

using RSLogic
using RCAdder
using RSGraphics
using LCommon

using Dice
using Plots
using Graphs

function get_model(parameters ::Dict{String, Any}) ::Dice.Model
    num_spins = length(parameters["problem"])
    if mod(num_spins - 1, 4) != 0
        error("Inconsistent problem")
    end
    
    num_bits = div(num_spins - 1, 4)
    return Dice.Model(RCAdder.get_rcadder_weighted_graph(num_bits),
                      Dice.model_2_hybrid_coupling,
                      parameters["delta_t"])
end

function get_onebit_model(parameters) ::Dice.Model
    return Dice.Model(RCAdder.generate_rcadder_graph(1),
                      Dice.model_2_hybrid_coupling,
                      parameters["delta_t"])
end

function get_agitation(num_spins ::Int) ::Vector{Float64}
    return RCAdder.get_random_x(num_spins)
end

function process_rca(model, parameters)
    displaced = parameters["displaced"]
    sig_res = length(displaced) + 1
    rights = zeros(Int, sig_res)
    for _ in 1:parameters["num_samples"]
        state, pinned, marks =
            RCAdder.get_initial_rca_rnd_marks(parameters["problem"], displaced)
        free = setdiff(1:nv(model.graph), pinned)

        for _ in 1:parameters["num_agitations"]
            state[2][free] = get_agitation(length(free))
            state = Dice.propagate_pinned(model,
                                          parameters["num_steps"],
                                          state,
                                          pinned)
        end
        
        stres = RCAdder.state_to_result_marks(state, displaced, marks)
        status = RCAdder.is_result_correct.(stres)
        rights .+= status[1:sig_res]
    end
    return rights
end

function get_slice_indices(slice ::Int) ::Vector{Int}
    stride = 5 # number of spins in FA
    return collect(1:stride) .+ (slice - 1) * (stride - 1)
end

function collapse_to_mark(x_val ::Float64, marks ::Vector{Float64})
    return marks[argmin(abs.(marks .- x_val))]
end

function collapse_small_state(small_state, marks)
    return (small_state[1],
            [collapse_to_mark(small_state[2][k], marks) for k in 1:length(small_state[2])])
end

function process_delayed_rca(model ::Dice.Model, parameters)
    log = 0
    displaced = parameters["displaced"]
    sig_res = length(displaced) + 1
    rights = zeros(Int, sig_res)
    for _ in 1:parameters["num_samples"]
        state, pinned, marks =
            RCAdder.get_initial_rca_rnd_marks(parameters["problem"], displaced)

        log > 0 && println("Marks: ", marks)

        for slice in 1:parameters["num_bits"]
            slice_inds = get_slice_indices(slice)
            small_state = (state[1][slice_inds], state[2][slice_inds])
            log > 0 && println("Initial: ", small_state)
            small_pinned = [1, 2, 3]
            small_state = Dice.propagate_pinned(model,
                                                parameters["num_steps"],
                                                small_state,
                                                small_pinned)
            log > 0 && println("Processed: ", small_state)
            # collapse small_state to the nearest mark
            small_state = collapse_small_state(small_state, marks)
            log > 0 && println("Collapsed: ", small_state)
            state[1][slice_inds] .= small_state[1]
            state[2][slice_inds] .= small_state[2]
        end
        
        stres = RCAdder.state_to_result_marks(state, displaced, marks)
        status = RCAdder.is_result_correct.(stres)
        rights .+= status[1:sig_res]

    end
    return rights
end


#########
##
## Main body
##
#########

x = 3411433493
y = 2079581652
cin = 0

fact = 41
parameters = Dict{String, Any}(
    "delta_t" => 0.1/fact,
    "total_time" => 2.0,
    "num_agitations" => 1,
    "num_samples" => 5000,
    "displaced" => Int[17, 3, 25, 0, 8, 30, 12, 21],
    "string_form" => "$x + $y + <$cin>",
)
parameters["num_steps"] = Int(round(parameters["total_time"] / parameters["delta_t"]))

prb_state = addition_to_spin_state(x, y, cin)
mod(length(prb_state) - 1, 4) != 0 &&
    error("Wrong addition problem")

parameters["problem"] = prb_state
parameters["num_bits"] = div(length(prb_state) - 1, 4)

model = get_model(parameters)

time_range = 6.1:0.65:31.7

out_probs ::Vector{Vector{Float64}} = []
for t in time_range
    println("Time: $t")
    parameters["total_time"] = t
    parameters["num_steps"] = Int(round(parameters["total_time"] / parameters["delta_t"]))
    rights = process_rca(model, parameters)
    println("Time: $t, rights: $rights out of $(parameters["num_samples"])")
    push!(out_probs, rights)
end

p = plot()

for k in 1:length(out_probs[1])
    dd = [out_probs[r][k]/parameters["num_samples"] for r in 1:length(out_probs)]
    plot!(time_range, dd, legend = false,
          marker = (:circle, 4),
          xlabel = "Machine time",
          ylabel = "Convergence probability",
          title = "x = $x, y = $y, cin = $cin", dpi = 200)
end

filename = "CONL_x-$(x)_y-$(y)_cin-$(cin)_s-$(fact)"
filename *= "_a-$(parameters["num_agitations"])"
filename *= "_i-$(length(parameters["displaced"]))"
filename *= "_f-$(fact))"
filename = sim_output_dir * filename

out_fig_file = filename * ".png"
savefig(out_fig_file)
println("Figure saved to $out_fig_file")

# Dump data

out_data_file = filename * ".dat"
open(out_data_file, "a") do out
    println(out, "# Convergence vs machine time (symmetric displacements)\n")
    println(out, "#PARAMETERS")
    for (key, value) in parameters
        println(out, "# $key : $value")
    end

    println(out, "\nseries : {")
    for (t, rights) in zip(time_range, out_probs)
        println(out, "\tpoint : {")
        println(out, "\t\ttime = $t")
        println(out, "\t\trights = $rights")
        println(out, "\t}")
    end    
    println(out, "}")
end
println("Data saved to $out_data_file")
