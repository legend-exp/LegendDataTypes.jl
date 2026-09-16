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
    @test_throws "Unknown physical unit \"not_a_unit\"" units_from_string("not_a_unit")
    @test_throws "Unknown physical unit" units_from_string("keV/")
    @test_throws "Unknown physical unit" units_from_string("keV + 1")
    @test units_from_string("ns^-1") == u"ns^-1"
    @test units_from_string("keV/s") == u"keV/s"
    @test units_from_string("m*s^-2") == u"m*s^-2"
    @test units_from_string("μs") == u"μs"
    for u in (u"keV", u"ns", u"mm", u"ns^-1", u"keV/s", NoUnits)
        @test units_from_string(units_to_string(u)) == u
    end

    # A unit string is parsed without `eval`, so it can be read in a precompile workload.
    mktempdir() do dir
        mkpath(joinpath(dir, "src"))
        write(joinpath(dir, "Project.toml"), """
            name = "UnitsInWorkload"
            uuid = "3a0f1b5e-8d51-4c8f-9a7e-5d6c2b1f0a11"
            version = "0.1.0"
            [deps]
            LegendDataTypes = "$(Base.PkgId(LegendDataTypes).uuid)"
            PrecompileTools = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
            """)
        write(joinpath(dir, "src", "UnitsInWorkload.jl"), """
            module UnitsInWorkload
            using LegendDataTypes, PrecompileTools
            @compile_workload LegendDataTypes.units_from_string("keV/s")
            end
            """)
        code = "import Pkg; Pkg.develop(path = \"$(pkgdir(LegendDataTypes))\"; io = devnull); " *
            "Pkg.add(\"PrecompileTools\"; io = devnull); using UnitsInWorkload"
        # The test sandbox leaves the stdlib out of the load path, where the subprocess needs Pkg.
        cmd = addenv(`$(Base.julia_cmd()) --startup-file=no --project=$dir -e $code`, "JULIA_LOAD_PATH" => "@:@stdlib")
        out = IOBuffer()
        ok = success(pipeline(cmd; stdout = out, stderr = out))
        ok || println(String(take!(out)))
        @test ok
    end
end
