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
            # for line in eachline(mtxFile) .|> strip .|> lowercase
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

Write the sparce adjacency matrix of `Graph` to file `filename` in the
reduced `MatrixMarket` format.

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

# Saving trajectories with metadata

using JSON3
using Dates

# const FVector = Vector{Float64}
# const SpinConf = Vector{Int8}
# const Hybrid = Tuple{SpinConf,FVector}

abstract type HeaderLine end

struct Comment <: HeaderLine
    text::String
end

struct Metadata{T} <: HeaderLine
    object::T
end

struct RawHeaderLine <: HeaderLine
    line::String
end


"""
    comment(text)

Create an ordinary human-readable header comment.

The resulting header entry is written as a line beginning with `# `.

# Example

```julia
comment("Simulation metadata")
```

writes

```text
# Simulation metadata
```
"""
comment(s::AbstractString) = Comment(String(s))


"""
    metadata(obj)

Create a machine-readable metadata header entry from an arbitrary
JSON-serializable object.

The object is written as a `#@` JSON record. Unlike [`record`](@ref), this
function does not add or require a `"type"` field.

# Example

```julia
metadata((git_commit = "abc123", seed = 42))
```

writes

```text
#@ {"git_commit":"abc123","seed":42}
```
"""
metadata(obj) = Metadata(obj)


"""
    raw_header_line(line)

Create a raw header line that is written exactly as supplied.

This is an escape hatch for callers who need full control over the emitted
header text. In normal use, prefer [`comment`](@ref), [`record`](@ref), or
[`metadata`](@ref).

# Example

```julia
raw_header_line("# already formatted")
```

writes

```text
# already formatted
```
"""
raw_header_line(s::AbstractString) = RawHeaderLine(String(s))


"""
    record(type; kwargs...)

Create a machine-readable metadata record with a `"type"` field.

The result is written as a `#@` JSON object. This is the recommended
constructor for both standard type-based metadata and user-defined
extension records.

# Examples

```julia
record(:node; index = 5, name = "pass_a")
```

writes

```text
#@ {"type":"node","index":5,"name":"pass_a"}
```

User-defined extension records are also allowed:

```julia
record(:solver; name = "my-relaxer", tolerance = 1e-8)
```
"""
function record(type::Union{Symbol,AbstractString}; kwargs...)
    if :type in keys(kwargs)
        throw(ArgumentError("do not pass `type` as a keyword; use the first argument to `record`"))
    end

    return Metadata((; type = String(type), kwargs...))
end


"""
    format_record(; name="relaxed-spins", version=1)

Create the standard relaxed-spin format declaration metadata record.

# Example

```julia
format_record()
```

writes

```text
#@ {"type":"format","name":"relaxed-spins","version":1}
```
"""
format_record(; name = "relaxed-spins", version = 1) =
    record(:format; name = String(name), version = Int(version))


"""
    created(value)
    created(dt::DateTime)

Create a standard `"created"` metadata record.

When given a string, `value` is written exactly as supplied. When given a
`DateTime`, the value is formatted as an ISO-like UTC timestamp ending in
`Z`.

# Examples

```julia
created("2026-01-26T18:44:00Z")
created(DateTime(2026, 1, 26, 18, 44, 0))
```
"""
created(value::AbstractString) =
    record(:created; value = String(value))


"""
    created(dt::DateTime)

Create a standard `"created"` metadata record from a `DateTime`.

The timestamp is formatted as

```text
yyyy-mm-ddTHH:MM:SSZ
```

This method assumes the supplied `DateTime` should be interpreted as UTC
and appends a trailing `Z`.
"""
function created(dt::DateTime)
    value = Dates.format(dt, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"
    return record(:created; value)
end

"""
    created()

Create a standard "created" metadata record using the current UTC time.
"""
function created()
    return created(Dates.now(Dates.UTC))
end

"""
    title(text)

Create a standard `"title"` metadata record.

# Example

```julia
title("Trajectory for test run 42")
```

writes

```text
#@ {"type":"title","value":"Trajectory for test run 42"}
```
"""
title(s::AbstractString) =
    record(:title; value = String(s))


"""
    notes(text)

Create a standard `"notes"` metadata record.

# Example

```julia
notes("Simulation after thermal relaxation.")
```

writes

```text
#@ {"type":"notes","value":"Simulation after thermal relaxation."}
```
"""
notes(s::AbstractString) =
    record(:notes; value = String(s))


"""
    node_indexing(base=1)

Create a standard `"node_indexing"` metadata record.

The `base` must be either `0` or `1`. Producers are encouraged to emit this
record explicitly to avoid ambiguity between indexing conventions.

# Examples

```julia
node_indexing(1)
node_indexing(0)
```
"""
function node_indexing(base::Integer = 1)
    base in (0, 1) ||
        throw(ArgumentError("node_indexing base must be 0 or 1"))

    return record(:node_indexing; base = Int(base))
end


"""
    node(index, name; kwargs...)

Create a standard `"node"` metadata record.

`index` is the logical node index, not a raw file column index. Extra keyword
arguments are included in the emitted JSON object, allowing user-defined
extension fields.

# Examples

```julia
node(1, "CTRL_1")
node(5, "pass_a"; role = "passivation")
```
"""
function node(index::Integer, name::AbstractString; kwargs...)
    return record(:node; index = Int(index), name = String(name), kwargs...)
end


"""
    as_header_line(x)

Normalize a user-supplied header item to an internal `HeaderLine`.

Supported inputs include:

- `HeaderLine` objects such as `comment(...)`, `record(...)`, and `metadata(...)`
- plain strings, treated as ordinary comments
- `NamedTuple`s, treated as machine-readable metadata
- `Dict`s and other `AbstractDict`s, treated as machine-readable metadata

This is primarily an internal helper used by [`write_header`](@ref).
"""
as_header_line(x::HeaderLine) = x


"""
    as_header_line(s::AbstractString)

Normalize a string header entry as an ordinary comment.
"""
as_header_line(s::AbstractString) =
    Comment(String(s))


"""
    as_header_line(nt::NamedTuple)

Normalize a `NamedTuple` header entry as machine-readable metadata.
"""
as_header_line(nt::NamedTuple) =
    Metadata(nt)


"""
    as_header_line(d::AbstractDict)

Normalize a dictionary header entry as machine-readable metadata.
"""
as_header_line(d::AbstractDict) =
    Metadata(d)


"""
    write_header_line(io, comment)

Write one ordinary comment header line to `io`.

If the comment text already begins with `#`, it is written unchanged. Otherwise,
`# ` is prepended.
"""
function write_header_line(io::IO, c::Comment)
    if startswith(c.text, "#")
        println(io, c.text)
    else
        println(io, "# ", c.text)
    end

    return nothing
end


"""
    write_header_line(io, metadata)

Write one machine-readable metadata header line to `io`.

The metadata object is serialized as JSON and prefixed with `#@ `.
"""
function write_header_line(io::IO, m::Metadata)
    print(io, "#@ ")
    JSON3.write(io, m.object)
    println(io)

    return nothing
end


"""
    write_header_line(io, raw)

Write one raw header line to `io` exactly as supplied.
"""
function write_header_line(io::IO, r::RawHeaderLine)
    println(io, r.line)
    return nothing
end


"""
    write_header(io, header)

Write all header entries to `io`.

Each item in `header` is normalized with [`as_header_line`](@ref) before being
written. This allows callers to mix explicit header constructors with plain
strings, `NamedTuple`s, and dictionaries.

# Example

```julia
header = [
    comment("Simulation metadata"),
    format_record(),
    node_indexing(1),
    node(1, "CTRL_1"),
    record(:solver; name = "my-relaxer"),
    metadata((git_commit = "abc123", seed = 42)),
]

write_header(stdout, header)
```
"""
function write_header(io::IO, header)
    for item in header
        write_header_line(io, as_header_line(item))
    end

    return nothing
end


"""
    default_time(traj)

Return the default time vector for `traj`.

The default time vector consists of 0-based integer indices:

```julia
0, 1, 2, ..., length(traj)-1
```

This helper is used by the `save_trajectories` methods that do not receive an
explicit `time` vector.
"""
default_time(traj::AbstractVector) =
    collect(0:(length(traj) - 1))


"""
    validate_trajectory_sample(spins, xs, k)

Validate one trajectory sample.

Checks that `spins` and `xs` have the same length and that every spin is either
`+1` or `-1`. The index `k` is used only for informative error messages.
"""
function validate_trajectory_sample(spins::AbstractVector, xs::AbstractVector, k)
    length(spins) == length(xs) ||
        throw(DimensionMismatch(
            "spin and FVector lengths differ at trajectory index $k: " *
            "$(length(spins)) != $(length(xs))"
        ))

    for s in spins
        if !(s == 1 || s == -1)
            throw(ArgumentError(
                "invalid spin value at trajectory index $k: $s; expected ±1"
            ))
        end
    end

    return nothing
end


"""
    save_trajectories(io::IO, traj::Vector{Hybrid}, time::AbstractVector{<:Real}, header=[])

Write a relaxed-spin trajectory file to `io`.

Each element of `traj` is a tuple

```julia
(spins, xs)::Tuple{Vector{Int8}, Vector{Float64}}
```

where `spins[i]` is the spin of node `i`, expected to be `+1` or `-1`, and
`xs[i]` is the corresponding floating-point value. The `time` vector must have
the same length as `traj`.

Each data line is written as

```text
time[k] spin_1 X_1 spin_2 X_2 ...
```

The optional `header` may contain ordinary comments and machine-readable
metadata. Supported header entries include:

- `comment("text")`, written as an ordinary `#` comment
- `record(:type; key=value, ...)`, written as a `#@` JSON metadata record
- `metadata(obj)`, where `obj` is any JSON-serializable object
- a `NamedTuple` or `Dict`, treated as machine-readable metadata
- a plain string, treated as an ordinary comment

# Example

```julia
traj = Hybrid[
    (Int8[1, -1, 1], Float64[0.12, 0.44, 0.91]),
    (Int8[1, -1, 1], Float64[0.15, 0.43, 0.90]),
]

time = [0.0, 0.1]

header = [
    comment("Simulation metadata"),
    format_record(),
    created("2026-01-26T18:44:00Z"),
    title("Trajectory for test run 42"),
    node_indexing(1),
    node(1, "CTRL_1"),
    node(3, "pass_a"),
    record(:solver; name = "my-relaxer", tolerance = 1e-8),
    metadata((git_commit = "abc123", seed = 42)),
]

save_trajectories(stdout, traj, time, header)
```
"""
function save_trajectories(
    io::IO,
    traj::Vector{Hybrid},
    time::AbstractVector{<:Real},
    header = HeaderLine[],
)
    length(traj) == length(time) ||
        throw(DimensionMismatch(
            "traj and time must have the same length: " *
            "$(length(traj)) != $(length(time))"
        ))

    write_header(io, header)

    expected_nodes = nothing

    for (k, ((spins, xs), t)) in enumerate(zip(traj, time))
        validate_trajectory_sample(spins, xs, k)

        if expected_nodes === nothing
            expected_nodes = length(spins)
        elseif length(spins) != expected_nodes
            throw(DimensionMismatch(
                "number of nodes changed at trajectory index $k: " *
                "$(length(spins)) != $expected_nodes"
            ))
        end

        print(io, t)

        for (s, x) in zip(spins, xs)
            print(io, ' ', Int(s), ' ', x)
        end

        println(io)
    end

    return nothing
end


"""
    save_trajectories(path::AbstractString, traj::Vector{Hybrid}, time::AbstractVector{<:Real}, header=[])

Write a relaxed-spin trajectory file to `path`.

This is a filename convenience wrapper around

```julia
save_trajectories(io, traj, time, header)
```

# Example

```julia
save_trajectories("trajectory.dat", traj, time, header)
```
"""
function save_trajectories(
    path::AbstractString,
    traj::Vector{Hybrid},
    time::AbstractVector{<:Real},
    header = HeaderLine[],
)
    open(path, "w") do io
        save_trajectories(io, traj, time, header)
    end

    return nothing
end


"""
    save_trajectories(traj::Vector{Hybrid}, time::AbstractVector{<:Real}, header=[])

Return a relaxed-spin trajectory file as a string.

This is a string-producing convenience wrapper around

```julia
save_trajectories(io, traj, time, header)
```

# Example

```julia
text = save_trajectories(traj, time, header)
```
"""
function save_trajectories(
    traj::Vector{Hybrid},
    time::AbstractVector{<:Real},
    header = HeaderLine[],
)
    return sprint() do io
        save_trajectories(io, traj, time, header)
    end
end


"""
    save_trajectories(io::IO, traj::Vector{Hybrid}, header=[])

Write a relaxed-spin trajectory file to `io` using default 0-based integer
times: : `0,  1, 2, ..., length(traj)-1`.

Thus the first data line begins with `0`, the second with `1`, and so on.

# Example

```julia
save_trajectories(stdout, traj, header)
```

is equivalent to

```julia
save_trajectories(stdout, traj, collect(0:(length(traj)-1)), header)
```
"""
function save_trajectories(
    io::IO,
    traj::Vector{Hybrid},
    header = HeaderLine[],
)
    time = default_time(traj)
    return save_trajectories(io, traj, time, header)
end


"""
    save_trajectories(path::AbstractString, traj::Vector{Hybrid}, header=[])

Write a relaxed-spin trajectory file to `path` using default 0-based integer
times: `0,  1, 2, ..., length(traj)-1`.

# Example

```julia
save_trajectories("trajectory.txt", traj, header)
```
"""
function save_trajectories(
    path::AbstractString,
    traj::Vector{Hybrid},
    header = HeaderLine[],
)
    time = default_time(traj)
    return save_trajectories(path, traj, time, header)
end


"""
    save_trajectories(traj::Vector{Hybrid}, header=[])

Return a relaxed-spin trajectory file as a string using default 0-based integer
times: `0,  1, 2, ..., length(traj)-1`

# Example

```julia
text = save_trajectories(traj, header)
```
"""
function save_trajectories(
    traj::Vector{Hybrid},
    header = HeaderLine[],
)
    time = default_time(traj)
    return save_trajectories(traj, time, header)
end


## Streaming versions

"""
    stream_header(io::IO, header)

Write metadata records to `io` as newline-delimited JSON.

This function is the stream-native counterpart of file header writing. It
emits machine-readable metadata records directly as JSON objects. If a
human-readable note should be transmitted through the stream, represent
it as metadata instead.

`header` may be either a collection of header entries or a single header
entry. Supported entries are the same as for file headers:

- `record(:type; key=value, ...)`
- `metadata(obj)`
- `NamedTuple`
- `Dict`

Ordinary comments are discarded. `RawHeaderLine` values are rejected.

# Example

```julia
header = [
    comment("This comment is not streamed"),
    format_record(),
    node_indexing(1),
    node(1, "CTRL_1"),
    record(:solver; name = "my-relaxer", tolerance = 1e-8),
]

stream_header(stdout, header)
```

writes

```text
{"type":"format","name":"relaxed-spins","version":1}
{"type":"node_indexing","base":1}
{"type":"node","index":1,"name":"CTRL_1"}
{"type":"solver","name":"my-relaxer","tolerance":1.0e-8}
```

The same function may be used later in a stream to interject updated
metadata:

```julia
stream_header(io, record(:conditions; temperature = 0.25, field = 1.5))
```
"""
function stream_header(io::IO, header)
    for item in header
        stream_header_item(io, item)
    end

    return nothing
end


"""
    stream_header(io::IO, item::HeaderLine)
    stream_header(io::IO, item::NamedTuple)
    stream_header(io::IO, item::AbstractDict)
    stream_header(io::IO, item::AbstractString)

Write a single metadata/header item to `io` as NDJSON.

This method is useful for interjecting one metadata record into an ongoing
stream.

# Example

```julia
stream_header(io, record(:conditions; coupling = 0.25))
```
"""
function stream_header(io::IO, item::HeaderLine)
    stream_header_item(io, item)
    return nothing
end

function stream_header(io::IO, item::NamedTuple)
    stream_header_item(io, item)
    return nothing
end

function stream_header(io::IO, item::AbstractDict)
    stream_header_item(io, item)
    return nothing
end

function stream_header(io::IO, item::AbstractString)
    stream_header_item(io, item)
    return nothing
end


"""
    stream_header(header) -> String

Return metadata records as a newline-delimited JSON string.

# Example

```julia
text = stream_header([
    format_record(),
    node_indexing(1),
    node(1, "CTRL_1"),
])
```
"""
function stream_header(header)
    return sprint() do io
        stream_header(io, header)
    end
end


"""
    stream_header_item(io::IO, item)

Write one normalized stream metadata item to `io`.

This is an internal helper used by [`stream_header`](@ref). Ordinary
comments are discarded. Raw header lines are rejected. Machine-readable
metadata is serialized as one JSON object followed by a newline.
"""
function stream_header_item(io::IO, item)
    h = as_header_line(item)

    if h isa Comment
        return nothing

    elseif h isa Metadata
        JSON3.write(io, h.object)
        println(io)
        return nothing

    elseif h isa RawHeaderLine
        throw(ArgumentError(
            "RawHeaderLine values cannot be streamed as NDJSON; " *
            "use metadata(...), record(...), a NamedTuple, or a Dict instead"
        ))

    else
        throw(ArgumentError("unsupported header item: $item"))
    end
end


"""
    sample_record(hybrid::Hybrid, time::Real; sample_index=nothing)

Construct the JSON-serializable stream record for a relaxed-spin configuration.

The emitted object has the form

```json
{"type":"sample","time":0.0,"r_spins":[{"state": [1,0.12]},{"state": [-1,0.44]}]}
```

where `r_spins[i] = [spin_i, X_i]`.

The sample is validated before the record is constructed. Each spin must be
`+1` or `-1`, and the spin and floating-point vectors must have equal length.

Currently, the optional `sample_index` is used only for informative error messages.
"""
function sample_record(
    hybrid::Hybrid,
    time::Real;
    sample_index = nothing,
)
    spins, xs = hybrid

    label = sample_index === nothing ? "streamed sample" : sample_index
    validate_trajectory_sample(spins, xs, label)

    r_spins = [(state = (Int(s), x),) for (s, x) in zip(spins, xs)]

    return (
        type = "sample",
        time = time,
        r_spins = r_spins,
    )
end


"""
    stream_sample(io::IO, hybrid::Hybrid, time::Real)

Write one relaxed-spin sample to `io` as a newline-delimited JSON object.

# Example

```julia
hybrid = (Int8[1, -1], Float64[0.12, 0.44])

stream_sample(stdout, hybrid, 0.0)
```

writes

```text
{"type":"sample","time":0.0,"r_spins":[{"state": [1,0.12]},{"state": [-1,0.44]}]}
```
"""
function stream_sample(io::IO, hybrid::Hybrid, time::Real)
    JSON3.write(io, sample_record(hybrid, time))
    println(io)

    return nothing
end


"""
    stream_sample(hybrid::Hybrid, time::Real) -> String

Return one relaxed-spin sample as a newline-delimited JSON string.

This is a string-producing convenience wrapper around `stream_sample`

# Example

```julia
text = stream_sample((Int8[1, -1], Float64[0.12, 0.44]), 0.0)
```
"""
function stream_sample(hybrid::Hybrid, time::Real)
    return sprint() do io
        stream_sample(io, hybrid, time)
    end
end


"""
    stream_trajectories(io::IO, traj::AbstractVector{<:Hybrid}, time::AbstractVector{<:Real}, header=[])

Write metadata and relaxed-spin trajectory samples to `io` as
newline-delimited JSON.

First, all JSON records from `header` are written using
[`stream_header`](@ref). Ordinary comments are discarded. Each
configuration in the trajectory is written as a `"sample"` message.

Each sample message has the form

```json
{"type":"sample","time":0.0,"r_spins":[{"state": [1,0.12]},{"state": [-1,0.44]}]}
```

The `traj` and `time` vectors must have the same length. The number of nodes
must remain constant across all samples. Each spin must be `+1` or `-1`.

# Example

```julia
traj = Hybrid[
    (Int8[1, -1], Float64[0.12, 0.44]),
    (Int8[1, -1], Float64[0.15, 0.43]),
]

time = [0.0, 0.1]

header = [
    comment("This comment is discarded in the stream"),
    format_record(),
    node_indexing(1),
    node(1, "CTRL_1"),
    node(3, "pass_a"),
]

stream_trajectories(stdout, traj, time, header)
```

writes

```text
{"type":"format","name":"relaxed-spins","version":1}
{"type":"node_indexing","base":1}
{"type":"node","index":1,"name":"CTRL_1"}
{"type":"node","index":3,"name":"pass_a"}
{"type":"sample","time":0.0,"r_spins":[{"state": [1,0.12]},{"state": [-1,0.44]}]}
{"type":"sample","time":0.1,"r_spins":[{"state": [1,0.15]},{"state": [-1,0.43]}]}
```

For live simulations, metadata can also be interjected manually between samples:
"""
function stream_trajectories(
    io::IO,
    traj::AbstractVector{<:Hybrid},
    time::AbstractVector{<:Real},
    header = HeaderLine[],
)
    length(traj) == length(time) ||
        throw(DimensionMismatch(
            "traj and time must have the same length: " *
            "$(length(traj)) != $(length(time))"
        ))

    stream_header(io, header)

    expected_nodes = nothing

    for (k, (hybrid, t)) in enumerate(zip(traj, time))
        spins, xs = hybrid

        record = sample_record(hybrid, t; sample_index = k)

        if expected_nodes === nothing
            expected_nodes = length(spins)
        elseif length(spins) != expected_nodes
            throw(DimensionMismatch(
                "number of nodes changed at trajectory index $k: " *
                "$(length(spins)) != $expected_nodes"
            ))
        end

        JSON3.write(io, record)
        println(io)
    end

    return nothing
end


"""
    stream_trajectories(io::IO, traj::AbstractVector{<:Hybrid}, header=[])

Write metadata and relaxed-spin trajectory samples to `io` as NDJSON, using
default 0-based integer times: `0, 1, 2, ..., length(traj)-1`.

# Example

```julia
stream_trajectories(stdout, traj, header)
```

is equivalent to

```julia
stream_trajectories(stdout, traj, collect(0:(length(traj)-1)), header)
```
"""
function stream_trajectories(
    io::IO,
    traj::AbstractVector{<:Hybrid},
    header = HeaderLine[],
)
    return stream_trajectories(io, traj, default_time(traj), header)
end


"""
    stream_trajectories(traj::AbstractVector{<:Hybrid}, time::AbstractVector{<:Real}, header=[]) -> String

Return metadata and relaxed-spin trajectory samples as a NDJSON string.

This is a string-producing convenience wrapper around `stream_trajectories`.

# Example

```julia
text = stream_trajectories(traj, time, header)
```
"""
function stream_trajectories(
    traj::AbstractVector{<:Hybrid},
    time::AbstractVector{<:Real},
    header = HeaderLine[],
)
    return sprint() do io
        stream_trajectories(io, traj, time, header)
    end
end


"""
    stream_trajectories(traj::AbstractVector{<:Hybrid}, header=[]) -> String

Return metadata and relaxed-spin trajectory samples as a newline-delimited JSON
string, using default 0-based integer times:: `0, 1, 2, ..., length(traj)-1`.

# Example

```julia
text = stream_trajectories(traj, header)
```
"""
function stream_trajectories(
    traj::AbstractVector{<:Hybrid},
    header = HeaderLine[],
)
    return sprint() do io
        stream_trajectories(io, traj, header)
    end
end

# One small usage pattern worth emphasizing

# stream_header(io, header)
# for step in 1:nsteps
#     hybrid, t = next_sample!(simulation)
#     stream_sample(io, hybrid, t)
#
#     if conditions_changed
#         stream_header(io, record(:conditions; coupling = T, field = B))
#     end
# end
