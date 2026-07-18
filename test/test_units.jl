# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

using Test
using LegendDataTypes
using Unitful, UnitfulAtomic

using LegendDataTypes: units_from_string, units_to_string

@testset "units from/to string" begin
    @test units_from_string("") == NoUnits
    @test units_from_string("none") == NoUnits
    @test units_from_string("ADC") == NoUnits
    @test units_from_string("adc") == NoUnits
    @test units_from_string("sample") == NoUnits
    @test units_from_string("keV") == u"keV"
    @test units_from_string("ns") == u"ns"
    @test units_from_string("e") == u"e_au"
    @test units_from_string("o/oo") == u"permille"
    @test units_from_string("o/o") == u"percent"
    @test_throws ArgumentError units_from_string("not_a_unit")
    for u in (u"keV", u"ns", u"mm", NoUnits)
        @test units_from_string(units_to_string(u)) == u
    end
end
