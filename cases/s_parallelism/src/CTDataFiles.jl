# module CTDataFiles
#
# functions for parsing data files containing the conversion rates vs
# machine time
#       Example file: TODO
#       Annotated Example: TODO

module CTDataFiles

export
    DataHolder,
    DHVector,
    processCTFile

const IVector = Vector{Int}
const FVector = Vector{Float64}
const SVector = Vector{String}

# Dice's types
const SpinConf = Vector{Int8}
const Hybrid = Tuple{SpinConf,FVector}

### Data types

"The data structure holding conversio vs time data. Corresponds to `point`
in the ct-data file"
mutable struct DataHolder
    time::Float64
    rights::IVector
end

DataHolder() = DataHolder(-1, [])
DataHolder(time::Float64, rights::Int64) = DataHolder(time, [rights])

"The vector of CT data points. Corresponds to `series` in the ct-data file"
const DHVector = Vector{DataHolder}

"The vector of CT data collections. Holds multiple `series` from the ct-data file"
const CTDataCollection = Vector{DHVector}

"The data type holding parameters"
const CTParameters = Dict{String, Any}

### String processing

"""
    string_split(str::AbstractString, sep = ' ') ::Tuple{String, String}

Split the string at the first occurrence of separator `sep` and return
`(car::String, cdr::String)` pair of the parts of the string strictly
before and strictly after `sep`. The returned parts are stripped of the
leading and trailing whitespaces. If the separator is at the beginning or
end of the non-white part of `str`, the respective returned parts are empty
strings. The default separator is space.
"""
function string_split(str::AbstractString, sep=' ')::Tuple{String,String}
    pos_sep = findfirst(sep, str)
    return pos_sep isa Nothing ?
           (str, "") :
           (strip(str[1:(pos_sep-1)]), strip(str[(pos_sep+1):end]))
end

"""
    split_string_at(str::AbstractString, samples::Vector{Char})

Split the string `str` at the first found separator from the list provided
in `samples` and return tuple `(before_separator, found_separator,
after_separator)` with `before_separator` and `after_separator` stripped.
If none of the separators are found, return `(str, "", "")`.
"""
function string_split_at(str::AbstractString,
    samples::Vector{Char})::Tuple{String,String,String}
    l_str = length(str)
    pos = l_str + 1
    for s in samples
        p = findfirst(s, str)
        p isa Nothing && continue

        pos = min(pos, p[1])
        break
    end
    if pos > l_str
        return (strip(str), "", "")
    end
    return (strip(str[1:(pos-1)]), str[pos:pos], strip(str[(pos+1):end]))
end

function string_split_at(str::AbstractString, sep ::Char)
    return string_split_at(str, [sep])
end

### FSM

abstract type State end

struct InvalidState <: State
    error_msg::String
end

struct Default <: State end

struct Series <: State end
struct finalizeSeries <: State end
struct Point <: State end
struct finalizePoint <: State end
struct Time <: State end
struct Rights <: State end

struct processControl <: State end
struct Parameter <: State end
struct finalizeParameters <: State end

struct Command <: State end
struct cmd_STOP <: State end

getElement(::CTDataCollection, ::State, ::AbstractString) =
    error("Undefined file parsing method")

function getElement(acc::CTDataCollection, state::Command, line::AbstractString)
    # process line starting with '#'
    # currently only END is implemented
    # TODO: parse the PARAMETERS block
    _, command = string_split(line, '#')
    command == "END" && return (cmd_STOP(), "")

    command == "PARAMETERS" && return (Parameter(), "")

    return (Default(), "")
end

function getElement(acc::CTDataCollection, state::Default, line::AbstractString)
    # assume a nonempty string without leading whitespace as input
    # and perform basic dispatching
    if line[1] == '#'
        return (Command(), line)
    end
    if line[1:6] == "series"  # currently equivalent to "series : {"
        push!(acc, DHVector())
        return (Series(), "")
    end

    return (InvalidState("Expected '#' or 'series', got '$line'"), "")
end

function getElement(acc::CTDataCollection, state::Series, line::AbstractString)
    # Expect
    #     1. "}" -> finalizeSeries
    #     2. "point : {" -> Point
    #     3. Invalid input
    line[1] == '}' && return (finalizeSeries(), line)

    data_type, rest = string_split(line, ':')

    if data_type == "point"
        push!(acc[end], DataHolder())
        return (Point(), "") # we drop the rest after :
    end

    return (InvalidState("Invalid element '$data_type' in line '$line' in the series group"), "")
end

function getElement(acc::CTDataCollection,
                    state::finalizeSeries,
                    line::AbstractString)
    return (Default(), "")
end

function getElement(acc::CTDataCollection, state::Point, line::AbstractString)
    # Expect
    #     1. "}" -> finalizePoint
    #     2. "rights = [" -> processRights
    #     3. "time = float" -> processTime
    #     4. Invalid input
    line[1] == '}' && return (finalizePoint(), line)

    data_type, rest = string_split(line, '=')

    data_type == "time" && return (Time(), rest)

    if data_type == "rights"
        # we consume '['
        _, rest = string_split(rest, '[')
        return (Rights(), rest)
    end

    return (InvalidState("Invalid element '$data_type' in line '$line' in the data point"), "")
end

function getElement(acc::CTDataCollection,
                    state::finalizePoint,
                    line::AbstractString)
    if acc[end][end].time == -1 || length(acc[end][end].rights) == 0
        return (InvalidState("Incomplete data point"), "")
    end
    return (Series(), "")
end

function getElement(acc::CTDataCollection, state::Time, line::AbstractString)
    # The positivity is the only constratint
    if acc[end][end].time != -1
        return (InvalidState("Multiple time specification in the data point"), "")
    end
    time = parse(Float64, line)
    if time <= 0
        return (InvalidState("Invalid nonpositive time is provided"), "")
    end
    acc[end][end].time = time
    return (Point(), "")
end

function getElement(acc::CTDataCollection,
                    state::Rights,
                    line::AbstractString)
     # Expect
    #     1. "]" -> check for nonempty and Point
    #     2. "," -> consume and processRights
    #     3. "int<,]>" -> check for positivity, push to the end processRights
    if line[1] == ']'
        length(acc[end][end].rights) == 0 && 
            return (InvalidState("Empty set of rights"), "")
        return (Point(), "")
    end
    line[1] == ',' && return (Rights(), line[2:end])
    car, sep, cdr = string_split_at(line, [',', ']'])
    val = parse(Float64, car)
    if val < 0
        return (InvalidState("Invalid negative value of rights counter"), "")
    end
    push!(acc[end][end].rights, val)
    return (Rights(), sep * cdr)
end

function getElement(parameters::CTDataCollection,
                    state::Parameter,
                    line::AbstractString)
    # lines without a content are taken as end of blocks
    if length(line) < 1 || line[1] != '#'
        return (finalizeParameters(), line)
    end
    # within blocks, lines with a content but without '=' or ":" considered
    # to be regular comments and ignored    
    _, sep, _ = string_split_at(line, ['=', ':'])
    isempty(sep) && return (Parameter(), "")
    _, _, line = string_split_at(line, ['#'])
    return (Parameter(), line)
end

# END of Control processing FSM

function processControls(block::SVector) ::Dict
    parameters = Dict()
    for line in block
        if !isempty(line)
            key, _, val = string_split_at(line, ['=', ':'])
            parameters[key] = val
        end
    end

    return parameters
end

### API functions

function processCTFile(file_name::String)::Tuple{CTDataCollection, Vector{Dict}}
    data_coll::CTDataCollection = []
    parameters_batches ::Vector{SVector} = []
    control_blocks ::SVector = []
    open(file_name, "r") do res_file
        pars_state = Default()
        for (linecount, line) in enumerate(eachline(res_file) .|> strip)
            while !isempty(line)
                (pars_state, line) = getElement(data_coll, pars_state, line)
                if pars_state isa Parameter && !isempty(line)
                    push!(control_blocks, line)
                    line = ""
                elseif pars_state isa finalizeParameters
                    push!(parameters_batches, control_blocks)
                    control_blocks = []
                    pars_state = Default()
                elseif pars_state isa InvalidState
                    @error "$file_name:$linecount: $(pars_state.error_msg)"
                end
            end
            pars_state isa cmd_STOP && break
        end
    end
    parameters = processControls.(parameters_batches)
    return (data_coll, parameters)
end

function processFile(file_name::String)::Tuple{DHVector, Vector{Dict}}
    return processCTFile(file_name::String)::Tuple{DHVector, Vector{Dict}}
end

end # module ends here
