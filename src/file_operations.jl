# Script file_operations.jl
#
# File operations

"""
    load_mtx_adjacency(filename::AbstractString,
                    min_weight::Float64 = weight_eps)::Matrix

Read the full adjacency matrix of the graph stored in file `filename`

This function reads the graph using `load_mtx_graph` and returns
its adjacency matrix converted to a full matrix.
"""
function load_mtx_adjacency(filename::AbstractString,
    min_weight::Float64=weight_eps)::Matrix
    return Matrix(adjacency_matrix(load_mtx_graph(filename, min_weight)))
end

"""
    load_mtx_graph(filename::AbstractString,
                min_weight::Float64 = weight_eps)::ModelGraph

Read graph stored in file `filename` in the simplified Matrixmarket
format. Weights with the magnitude below `min_weight` are taken to be 0.

# Output
    SimpleGraph or SimpleWeightedGraph object

# The expected file format
The simplified format:
    Header: lines starting with % (such lines are ignored)
    Descriptor:
             The valid Matrix Market descriptors are
             R C E - number of rows, columns, non-zero entries
             N E - number of nodes, edges

    Content
        u v w_{u,v}

## Example of an MTX file
% K_3 in the mtx format
3 3
1 2
1 3
2 3
"""
function load_mtx_graph(filename::AbstractString,
    min_weight::Float64=weight_eps)::ModelGraph
    open(filename, "r") do mtxFile
        header = true
        linecount = 0
        shape = [0, 0]
        while !eof(mtxFile)
            #            for line in eachline(mtxFile) .|> strip .|> lowercase
            line = readline(mtxFile) .|> strip .|> lowercase
            linecount += 1
            if isempty(line) || line[1] == '%'
                continue
            end
            tokens = split(line)
            # the first content holding line is the shape descriptor
            if length(tokens) < 2 || length(tokens) > 3
                error("In $filename:$linecount, the descriptor $line is invalid")
            end
            shape .+= parse.(Int64, tokens[[1, end]])
            header = false
            break
        end
        if header
            error("In $filename, no header was found")
        end

        sources = zeros(Int32, shape[2])
        dests = zeros(Int32, shape[2])
        weights = zeros(Float64, shape[2])
        idx = 1

        while !eof(mtxFile)
            #            for line in eachline(mtxFile) .|> strip .|> lowercase
            line = readline(mtxFile) .|> strip .|> lowercase
            linecount += 1
            if isempty(line) || line[1] == '%'
                continue
            end
            tokens = split(line)

            if length(tokens) < 2 || length(tokens) > 3
                error("In $filename:$linecount, content line $line is invalid")
            end
            u, v = parse.(Int64, tokens[1:2])
            w = length(tokens) == 3 ? parse(Float64, tokens[3]) : 1.0
            # we skip edges with too small weight
            abs(w) < min_weight && continue

            sources[idx] = u
            dests[idx] = v
            weights[idx] = w
            idx += 1
        end
        if maximum(weights) - minimum(weights) > weight_eps
            # we have a weighted graph
            return SimpleWeightedGraph(sources, dests, weights)
        end
        G = SimpleGraph(shape[1])
        for (u, v) in zip(sources, dests)
            Graphs.add_edge!(G, u, v)
        end
        return G
    end
end

"""
    loadDumpedGraph(filename::AbstractString)::SimpleGraph

NOTE: legacy alias for `loadMTXGraph`
"""
function loadDumpedGraph(filename::AbstractString)::ModelGraph
    @warn "Function `loadDumpedGraph` is obsolete use `load_mtx_graph` instead"
    return load_mtx_graph(filename)
end

"""
    save_mtx_graph(Graph::SimpleGraph, filename::String)

Writes the sparce adjacency matrix of `Graph` to file `filename` in the reduced
`MatrixMarket` format.

The adjacency matrix is written in the reduced `MatrixMarket` format.
The first line contains two numbers: `|V| |E|`
Next lines contain three numbers: `u v A_{u, v}`
where `u` and `v` are the nodes numbers and `A_{u, v}` is the edge weight.
"""
function save_mtx_graph(G::ModelGraph, filename::AbstractString)
    open(filename, "w") do out
        println(out, "$(nv(G)) $(ne(G))")
        for edge in edges(G)
            if G isa SimpleWeightedGraph
                w = weight(edge)
            else
                w = 1
            end
            println(out, "$(edge.src) $(edge.dst) $w")
        end
    end
end

function dumpGraph(G::ModelGraph, filename::AbstractString)
    @warn "Function `dumpGraph` is obsolete use `save_mtx_graph` instead"
    save_mtx_graph(G, filename)
end

"""
    dump_graph_console(G::ModelGraph)

Prints the sparce adjacency matrix of `Graph` to console
in the reduced `MatrixMarket` format (see `save_mtx_graph`).
"""
function dump_graph_console(G::ModelGraph)
    println("$(nv(G)) $(ne(G))")
    for edge in edges(G)
        if G isa SimpleWeightedGraph
            w = weight(edge)
        else
            w = 1
        end

        println("$(edge.src) $(edge.dst) $w")
    end
end

function dump_distribution(P, filename)
    # P is a tuple (x, Px) as, for instance, produced by the integer distribution
    open(filename, "w") do out
        for i in 1:length(P[1])
            println(out, P[1][i], " ", P[2][i])
        end
    end
end
