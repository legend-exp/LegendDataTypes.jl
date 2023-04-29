# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

"""
    TableLike = AbstractVector{<:NamedTuple}

The type must also support the Tables.jl interface.

Examples are `TypedTables.Table` and `StructArrays.StructArray`
"""
const TableLike = AbstractVector{<:NamedTuple}


@enum DAQType daq_gerda=0x38bce998
export DAQType, daq_gerda

@enum EventType evt_undef=0 evt_real=1 evt_pulser=2 evt_mc=3 evt_baseline=4
export EventType, evt_undef, evt_real, evt_pulser, evt_mc, evt_baseline


const daq_noflags = UInt32(0)
const daq_underflow = UInt32(1<<0)
const daq_overflow = UInt32(1<<1)
const daq_corrupt = UInt32(1<<2)
const daq_reserved1 = UInt32(1<<3)
const daq_prepileup = UInt32(1<<4)
const daq_postpileup = UInt32(1<<5)

export daq_noflags, daq_underflow, daq_overflow, daq_corrupt, daq_reserved1, daq_prepileup, daq_postpileup
