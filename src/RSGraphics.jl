# project: logic
# module RSGraphics.jl
#
# General functions for plotting relaxed spin configurations

module RSGraphics

using Plots
using Graphs
using SimpleWeightedGraphs

using LaTeXStrings

export
    make_circle,
    add_state,
    show_state,
    make_dash  # don't know why it doesn't work without it

const FMatrix = Matrix{Float64}
const IVector = Vector{Int}
const I8Vector = Vector{Int8}
const FVector = Vector{Float64}
const SVector = Vector{String}

### Graphics primitives

ang = range(0, 2π, length=60)
circ_radius = 1.0
arrow_size = 0.5

circle(x, y, r = circ_radius) = Plots.Shape(r * sin.(ang) .+ x, r * cos.(ang) .+ y)

arrow_up(x, y, ar_len = arrow_size) =
    Plots.plot!([x, x], [y - ar_len / 2, y + ar_len / 2],
                arrow=true, color=:blue, linewidth=3, label="")
arrow_down(x, y, ar_len = arrow_size) =
    Plots.plot!([x, x], [y + ar_len / 2, y - ar_len / 2],
                arrow=true, color=:blue, linewidth=3, label="")

arrow_down(x::Tuple{Float64, Float64}) = arrow_down(x[1], x[2])

get_pos(X ::Float64, radius::Float64 = circ_radius) =
    radius .* (sin(π * X), -cos(π * X))

set_arrow(X ::Float64, s ::Int8, ar_len::Float64 = 0.5) = begin
	ar_size = ar_len
	(x_pos, y_pos) = get_pos(X)
	Plots.plot!([x, x], [y - s * ar_size/2, y + s * ar_size/2], arrow = true, color = :blue, linewidth = 3, label = "")
end

# plots a vertical dash at point with the circle coordinate `X`.
# Only useful to show the zero point or the spin reversion boundary
show_dash(X ::Float64, dash_len ::Float64, c = :black) = begin
    (x_pos, y_pos) = get_pos(X)
    Plots.plot!([x_pos, x_pos],
                [y_pos - dash_len/2, y_pos + dash_len/2],
                color = c, linewidth = 2, label = false)
end

spin_point(X ::Float64) = begin 
	rad = 0.1
	(x_pos, y_pos) = get_pos(X)
	Plots.plot!(circle(x_pos, y_pos, rad), fillcolor = :grey, label = false)
end

show_gen_spin(X ::Float64, s ::Int8, ar_len ::Float64, color ::Symbol) = begin
    ar_size = ar_len
    rad = 0.09
    (x_pos, y_pos) = get_pos(X)
    Plots.plot!(circle(x_pos, y_pos, rad), fillcolor = :grey, label = false)
    Plots.plot!([x_pos, x_pos],
                [y_pos - 0 * s * ar_size/2, y_pos + s * ar_size/2],
                arrow = true, color = color, linewidth = 3, label = false)
end

show_main_spin(X ::Float64, s ::Int8) = show_gen_spin(X, s, 0.5, :blue)

show_central_spin(X ::Float64, s ::Int8) = show_gen_spin(X, s, 0.7, :red)

show_spin(X ::Float64, s ::Int8, selector ::Symbol) =
    selector == :central ? show_central_spin(X, s) : show_main_spin(X, s)

make_circle(radius = circ_radius, cut_off = 1.4) =
    Plots.plot(circle(0, 0, radius),
               xlim=(-cut_off, cut_off),
               ylim=(-cut_off, cut_off),
               fillcolor=:white,
               linewidth=2,
               label=false,
               ratio=1,
               grid = false,
               axis = ([], false))

function add_state(state ::Tuple{I8Vector, FVector})
    return add_state(state[1], state[2])
end

function add_state(Ss ::I8Vector, Xs ::FVector)
    num_spins = length(Ss)
    for i in 1:num_spins
        show_spin(Xs[i], Ss[i], i == num_spins ? :main : :main)
    end
end

function show_state(Ss ::I8Vector, Xs ::FVector)
    p1 = make_circle()
    add_state(Ss, Xs)
    show_dash(0.0, 0.1)
    return p1
end

function show_state(state ::Tuple{I8Vector, FVector})
    return show_state(state[1], state[2])
end

## Example

# ss ::Vector{Int8} = [1, 1, -1, -1]
# xs ::Vector{Float64} = [-0.2, 0.7, 0.34, -0.5]
#
# p1 = RSGraphics.make_circle()
# RSGraphics.add_state(ss, xs)
# RSGraphics.show_dash(0.0, 0.1)
# The last command is needed (for unknown reasons) to update the plot
#
# Or there's a function doing all of the above
# p1 = show_state(ss, xs)

# To show multiple states side by side, we can use the returned plots
# plot(p1, p2, p3, layout = (1,3), legend = false)

# To save to a png-file, it is useful to specify dpi, otherwise
# the figure may be too small
# plot(p1, p01, p2, layout = (1,3), legend = false, dpi = 200)
# file_out = "k$selector-spins.png"
# savefig(file_out)
# println("The figure is saved to '$file_out'")

end # module ends here
