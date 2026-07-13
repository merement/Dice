# Collection of obsolete functions
# Kept for historical reasons to revive some old scripts (not guaranteed)

# Source: src/statistical_methods.jl
function c_variance(V, intervals)
    # OBSOLETE (see NOTE)
    # Calculates first two momenta of the pieces of data that fall
    # within given M intervals provided in `intervals`
    #
    # INPUT:
    #   V - array of data
    #   intervals - M x 2 array of intervals boundaries
    #
    # OUTPUT:
    #   M x 3 array with the number of points, mean and variance of data
    #   inside the respective intervals
    #
    # NOTE: cannot process a single interval and does not treat folding intervals

    out = zeros(size(intervals)[1], 3)

    for i in 1:(size(out)[1])
        ari = filter(t -> intervals[i,1] < t < intervals[i,2], V)
        ni = length(ari)
        av = sum(ari) / ni
        var = sum((ari .- av).^2) / ni

        out[i, :] = [ni, av, sqrt(var)]
    end

    return out
end
