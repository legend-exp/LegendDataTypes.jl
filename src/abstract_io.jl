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
    readdata(input, name::AbstractString)
    readdata(input, name::AbstractString, SomeDataType::Type)

Read the value stored under `name` from `input`.

LEGEND I/O packages add methods for the I/O-object types they handle.
"""
function readdata end


"""
    writedata(output, name::AbstractString, x)

Write a value `x` under `name` to `output`.

LEGEND I/O packages add methods for the I/O-object types they handle.
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

@inline getunits(::Nothing) = NoUnits


"""
    setunits!(x)

Set the units of x.

LEGEND I/O packages will need to add methods for the I/O-object types they
handle.
"""
function setunits! end

const _pseudo_unit_expr = r"(^|[^A-Za-z])(ADC|adc|sample)([^A-Za-z]|$)"

# Unit strings that name something else in the unit modules, or nothing at all.
const _unit_aliases = Dict("e" => u"e_au", "o/oo" => u"permille", "o/o" => u"percent")

# A unit expression is evaluated here rather than with `eval`, which cannot run while a
# package reading units in its precompile workload is precompiled. `Unitful.lookup_units`
# has replaced the unit symbols by their units and admits only arithmetic on them.
_eval_units(x) = x
_eval_units(ex::Expr) = getfield(Base, ex.args[1])(map(_eval_units, ex.args[2:end])...)

function units_from_string(s::AbstractString)
    if isempty(s) || s == "none"
        NoUnits
    elseif !isnothing(match(_pseudo_unit_expr, s))
        NoUnits
    elseif haskey(_unit_aliases, s)
        _unit_aliases[s]
    else
        try
            _eval_units(Unitful.lookup_units([Unitful, UnitfulAtomic], Meta.parse(s)))
        catch
            throw(ArgumentError("Unknown physical unit \"$s\""))
        end
    end
end

units_to_string(u::Unitful.Unitlike) = string(u)
