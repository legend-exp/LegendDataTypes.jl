# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

module LegendDataTypes

using Printf

using ArraysOfArrays
using ElasticArrays
using EncodedArrays
using IntervalSets
using RecipesBase
using StaticArrays
using StatsBase
using StructArrays
using Unitful
using UnitfulAtomic
using UnitfulParsableString

import Tables
import TypedTables
import Colors

using RadiationDetectorSignals: RDWaveform

include("daq.jl")
include("type_registry.jl")
include("array_codecs.jl")
include("abstract_io.jl")
include("output_generation.jl")
include("radware_sigcompress.jl")
include("utils.jl")
include("data_filters.jl")


const array_codecs = TypeRegistry{AbstractArrayCodec}()


function __init__()
    array_codecs[:uleb128_zigzag_diff] = VarlenDiffArrayCodec
    array_codecs[:radware_sigcompress] = RadwareSigcompress

    # ToDo: Remove, deprecated:
    array_codecs[:varlend_diff1] = VarlenDiffArrayCodec
end

end # module
