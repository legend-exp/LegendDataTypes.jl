# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).


# supports getindex (and possibly view) and close
abstract type AbstractLegendInput end
export AbstractLegendInput

# support setindex!, append! and close
abstract type AbstractLegendOutput end
export AbstractLegendOutput



struct LegendNullInput <: AbstractLegendInput end
export LegendNullInput


function Base.open(filename::AbstractString, ::Type{LegendNullInput})
    LegendNullInput()   
end

Base.close(input::LegendNullInput) = nothing

"""
    getindex(input::AbstractLegendInput, key::AbstractString)
    getindex(input::AbstractLegendInput, key::AbstractString, idxs::AbstractVector)

Get object at `key` from input.
"""
function Base.getindex(input::AbstractLegendInput, key::AbstractString, args...)
    throw(ErrorException("No key \"$key\" found in input $input"))
end

# TODO: Base.view



struct LegendNullOutput <: AbstractLegendOutput end
export LegendNullOutput

function Base.open(filename::AbstractString, ::Type{LegendNullOutput})
    LegendNullOutput()
end

Base.close(input::LegendNullOutput) = nothing

"""
    setindex!(output::AbstractLegendOutput, key::AbstractString)
    getindex(output::AbstractLegendOutput, key::AbstractString, idxs::AbstractVector)

Get object at `key` from input.
"""
function Base.getindex(output::LegendNullOutput, value::Any, key::AbstractString, args...)
    nothing
end

# TODO: Base.append!



# TODO: move following to the individual IO packages?


"""
    readdata(input, SomeDataType::Type)

Read a value of type `SomeDataType` from `input`.
"""
function readdata end


"""
    readdata(input, x::SomeDataType)

Write a value `x` to `output`.
"""
function writedata end


"""
    getunits(x)

Get the units of x, falls back to `Unitful.unit(x)` if no specialized method is
defined for the type of `x`.

LEGEND I/O packages shoud add methods for the I/O-object types they handle.
"""
function getunits end

@inline getunits(x::Any) = Unitful.unit(x)

Base.@pure getunits(::Nothing) = NoUnits


"""
    setunits!(x)

Set the units of x.

LEGEND I/O packages will need to add methods for the I/O-object types they
handle.
"""
function setunits! end


function units_from_string(s::AbstractString)
    if isempty(s)
        NoUnits
    else
        try
            Unitful.lookup_units(Unitful, Meta.parse(s))
        catch e
            if e isa ErrorException
                rethrow(ArgumentError("Unknown physical unit \"$s\""))
            else
                rethrow(e)
            end
        end
    end
end

units_to_string(u::Unitful.Unitlike) = "$u"
