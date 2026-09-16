# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

# The flatten of the column types the LEGEND data tiers hold, compiled once here: each type
# a table holds otherwise costs a specialization on first use, and a tier holds about ten.
@compile_workload begin
    for T in (Bool, Float32, Float64, Int16, Int32, Int64, UInt16, UInt32,
              typeof(1.0u"s"), typeof(1.0u"ns"), typeof(1.0u"μs"), typeof(1.0u"keV"), typeof(1.0u"ns^-1"), typeof(1.0u"e_au"),
              typeof(1.0f0u"ns"), typeof(1.0f0u"μs"), typeof(1u"ns"))
        v = [T[zero(T), zero(T)], T[zero(T)]]
        fast_flatten(v)
        flatten_by_key([(a = first(v), b = first(v)), (a = last(v), b = last(v))])
        sa = StructArray(a = first(v), b = first(v))
        fast_flatten([sa, sa])
        fast_flatten([TypedTables.Table(sa), TypedTables.Table(sa)])
    end
    foreach(units_from_string, ("", "none", "ADC", "keV", "ns", "μs", "s", "ns^-1", "e", "o/o"))
end
