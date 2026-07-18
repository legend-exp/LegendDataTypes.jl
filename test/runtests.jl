# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

import Test
Test.@testset "Package LegendDataTypes" begin

include("test_units.jl")
include("test_radware_sigcompress.jl")

end # testset
