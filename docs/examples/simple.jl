# script simple.jl
# Simple IM demonstration

using Dice

using Graphs
using SimpleWeightedGraphs

# Create a graph
graph = Graphs.Graph(7)
add_edge!(graph, 1, 2)
add_edge!(graph, 2, 3)
add_edge!(graph, 3, 4)
add_edge!(graph, 4, 5)
add_edge!(graph, 5, 6)
add_edge!(graph, 6, 7)
add_edge!(graph, 7, 1)

# Create a model
# Stop time in model units
total_time = 4.4
# Number of timesteps
num_steps = 100
delta_t = total_time/num_steps
model = Dice.Model(graph, Dice.model_2_hybrid_coupling, delta_t)

num_vertices = Graphs.nv(model.graph)
# A specific initial state
state = (Int8[-1, 1, 1, -1, -1, 1, 1], [0.117, 0.557, -0.164, 0.122, 0.471, -0.323, 0.224])
# Alternatively, one can generate a random initial state
# state = Dice.get_random_hybrid(num_vertices, 1.2)

println("Initial state")
@show state
state = Dice.propagate(model, num_steps, state)
println("\nFinal binary component of the R-spins: $(state[1])") # should be Int8[1, 1, -1, 1, -1, 1, -1]
cut = Dice.cut(model, state)
println("Final Cut value: $cut") # should be 6
