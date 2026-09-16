# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

using Test
using LegendDataTypes
using StructArrays, TypedTables, Unitful

using LegendDataTypes: fast_flatten, flatten_by_key

@testset "fast_flatten" begin
    # a wide table of few column types, flattened column by column
    names = Tuple(Symbol("c$i") for i in 1:60)
    cols(n) = NamedTuple{names}(Tuple(i % 3 == 0 ? rand(Int32, n) : i % 3 == 1 ? rand(n) : rand(n) .* u"keV" for i in 1:60))
    parts = [cols(3), cols(0), cols(5)]
    nt = flatten_by_key(parts)
    @test nt isa NamedTuple{names}
    @test all(nt[k] == vcat((p[k] for p in parts)...) for k in names)
    @test all(typeof(nt[k]) == typeof(parts[1][k]) for k in names)

    sas = StructArray.(parts)
    sa = fast_flatten(sas)
    @test sa isa StructArray && eltype(sa) == eltype(first(sas))
    @test sa == vcat(sas...)
    @test length(sa) == 8

    tbls = Table.(parts)
    tbl = fast_flatten(tbls)
    @test tbl isa Table && tbl == vcat(tbls...)

    # rows of a single part come back unchanged
    @test fast_flatten([sas[1]]) == sas[1]
end
