# Statistical methods

"""
    integer_distribution(arr::AbstractVector{Int})::Tuple{Vector{Int}, Vector{Float64}}

Evaluate a discrete (over integers) distribution function inferred from
the provided sample `arr`.

For each integer in [minimum(arr), maximum(arr)], the function counts
the relative number of incidences.

# INPUT:
    `arr::AbstractVector{Int}` - empirical sample of the distribution function

# OUTPUT:
    `(x::Vector{Int}, p(x)::Vector{Float64})`, where
    `x` - values sampled in `arr`
    `p(x)` - occurrence probabilities
"""
function integer_distribution(arr::AbstractVector{Int})::Tuple{Vector{Int}, Vector{Float64}}
    sorted_arr = sort(arr)
    n = length(sorted_arr)
    if n == 0
        return (Int[], Float64[])
    end

    vals = Int[]
    pvals = Float64[]
    
    sizehint!(vals, n)
    sizehint!(pvals, n)

    i = 1
    while i <= n
        val = sorted_arr[i]
        count = 1
        while i + count <= n && sorted_arr[i + count] == val
            count += 1
        end
        push!(vals, val)
        push!(pvals, count / n)
        i += count
    end

    return (vals, pvals)
end

"""
    cumulative_integer_distribution(arr::AbstractVector{Int})::Tuple{Vector{Int}, Vector{Float64}}

Evaluate a cumulative discrete distribution function over integers inferred
from the samples provided in `arr`.

For each integer `val` in [minimum(arr), maximum(arr)], the function counts
the relative number of values in `arr` strictly smaller than `val`.

# INPUT:
    `arr::AbstractVector{Int}` - empirical sample of the distribution function

# OUTPUT:
    `(x::Vector{Int}, p(x)::Vector{Float64})`, where
    ``x`` - values sampled in `arr`
    ``p(x) = P(a < x)`` - probability to encounter value smaller than `x`
"""
function cumulative_integer_distribution(arr::AbstractVector{Int})::Tuple{Vector{Int}, Vector{Float64}}
    sorted_arr = sort(arr)
    n = length(sorted_arr)
    if n == 0
        return (Int[], Float64[])
    end

    vals = Int[sorted_arr[1] - 1]
    pvals = Float64[0.0]
    
    sizehint!(vals, n + 2)
    sizehint!(pvals, n + 2)

    i = 1
    while i < n
        val = sorted_arr[i]
        
        j = i + 1
        while j <= n && sorted_arr[j] == val
            j += 1
        end
        
        if j > n
            push!(vals, val + 1)
            push!(pvals, Float64(n))
            return (vals, pvals ./ n)
        end
        
        push!(vals, val)
        push!(pvals, pvals[end] + (j - i))
        i = j
    end

    push!(vals, sorted_arr[i] + 1)
    push!(pvals, Float64(n))
    return (vals, pvals ./ n)
end

"""
    av_dispersion(graph::SimpleGraph, V::AbstractVector{<:Real})::Float64

Evaluates average difference between the dynamical variables at adjacent nodes.

# INPUT:
    `graph::SimpleGraph`
    `V::AbstractVector{<:Real}` - the array of dynamical variables

# OUTPUT:
    `Σ_e |V(a) - V(b)| / |graph|`
"""
function av_dispersion(graph::SimpleGraph, V::AbstractVector{<:Real})::Float64
    return sum(abs(V[edge.src] - V[edge.dst]) for edge in edges(graph)) / nv(graph)
end

"""
    cluster_variance(V::AbstractVector{<:Real}, interval)::Vector{Float64}

Calculate first two momenta of the pieces of data that fall within
the provided interval.
    
# INPUT:
    `V::AbstractVector{<:Real}` - array of data
    `interval` - (vmin, vmax)
    
# OUTPUT:
    `Array[3] = [number of points, mean, variance]`
"""
function cluster_variance(V::AbstractVector{<:Real}, interval::Union{Tuple{<:Real, <:Real}, AbstractVector{<:Real}})::Vector{Float64}
    vmin, vmax = interval[1], interval[2]
    
    count = 0
    sum_val = 0.0
    for t in V
        if vmin < t < vmax
            count += 1
            sum_val += t
        end
    end
    
    if count == 0
        return Float64[0.0, 0.0, 0.0]
    end
    
    av = sum_val / count
    
    sum_sq_diff = 0.0
    for t in V
        if vmin < t < vmax
            sum_sq_diff += (t - av)^2
        end
    end
    
    var = sum_sq_diff / count
    return Float64[count, av, sqrt(var)]
end
