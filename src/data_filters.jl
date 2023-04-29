# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).


const _ch_string_exp = r"ch([0-9]+)$"

"""
    chname2int(channel_string::AbstractString)::Integer

Convert a channel name string, as used in LEGEND data files to an integer
channel number.
"""
function chname2int(channel_string::AbstractString)
    parse(Int, match(_ch_string_exp, channel_string)[1])
end
export chname2int


"""
    int2chname(channel_number::Integer)::AbstractString

Convert an integer channel number to a channel name string as used in LEGEND
data files.
"""
int2chname(ch_no) = @sprintf("ch%03d", ch_no)
export int2chname


"""
    get_all_channels(ds::AbstractDict{<:AbstractString,<:Any})::AbstractVector{<:Integer}

Get the channel numbers for all channels in `datastore`.

Channels are identified as keys starting with "ch...", according to LEGEND
convention. `datastore` will typically be on-disk, e.g. a
`LegendHDF5IO.LegendHDFLHDataStore`.

Channel named are mapped to integer channel numbers via [`chname2int`](@ref).
"""
function get_all_channels(datastore::AbstractDict{<:AbstractString,<:Any})
    sort(chname2int.(filter(startswith("ch"), keys(datastore))))
end
export get_all_channels


"""
    get_raw_ch_data(ds::AbstractDict{<:AbstractString,<:Any}, ch::Integer)::TableLike

Get the raw data table for channel `ch` in `datastore`.

`datastore` will typically be on-disk, e.g. a
`LegendHDF5IO.LegendHDFLHDataStore`. Channel numbers are mapped to channel
names via [`int2chname`](@ref).
"""
function get_raw_ch_data(ds::AbstractDict{<:AbstractString,<:Any}, ch::Integer)
    (ds[int2chname(ch)].raw)::TableLike
end
export get_raw_ch_data


"""
    get_daqenergy(datastore::AbstractDict{<:AbstractString,<:Any}, ch::Integer)::AbstractVector{<:Integer}

Get the daq energy reconstruction contained in the raw data of channel `ch`
in `datastore`.

`datastore` will typically be on-disk, e.g. a `LegendHDF5IO.LegendHDFLHDataStore`.
Channel numbers are mapped to channel names via [`int2chname`](@ref).
"""
function get_daqenergy(datastore::AbstractDict{<:AbstractString,<:Any}, ch::Integer)
    if hasfield(typeof(datastore), :data_store)
        # Default implementation is slower than it could be for
        # LegendHDF5IO.LHDataStore, so use workaround.
        # ToDo: Fix in LegendHDF5IO.
        datastore.data_store[int2chname(ch)]["raw"]["daqenergy"][:]
    else
        get_raw_ch_data(datastore, ch).daqenergy[:]
    end::AbstractVector{<:Integer}
end
export get_daqenergy



export filter_raw_data_by_energy

"""
    filter_raw_data_by_energy(
        raw_data::TableLike,
        calib_func::Function,
        energy_windows::AbstractDict{Symbol,<:AbstractInterval{<:Number}}
    )::IDDict{Symbol,<:Any}

Filter the table `raw data`, selecing only events in the energy intervals in
`values(energy_windows)`.

The selection is based on DAQ energy reconstruction and the energy
calibration function `calib_func`.

Returns a dicts of raw data tables with the same keys as `energy_windows`.
"""
function filter_raw_data_by_energy(
    raw_data::TableLike,
    calib_func::Function,
    energy_windows::AbstractDict{Symbol,<:AbstractInterval{<:Number}};
    chunk_size = 10000
)
    chunk_idxs = Iterators.partition(eachindex(raw_data), chunk_size)
    flatten_by_key([let data = raw_data[idxs]
        cal_energy = calib_func.(data.daqenergy)
        IdDict((
            window_label => data[findall(x -> x in energy_windows[window_label], cal_energy)]
            for window_label in keys(energy_windows)
        ))
    end for idxs in chunk_idxs])
end


#=
"""
    filter_raw_data_by_energy(
        datastore::AbstractDict{<:AbstractString,<:Any},
        calib_funcs::IdDict{<:Integer,<:Function},
        energy_windows::AbstractDict{Symbol,<:AbstractInterval{<:Number}},
        chunk_size = 10000
    )::IDDict{Int,<:IDDict{Symbol,<:Any}}

Apply `filter_raw_data_by_energy` for multiple channels.

Uses the keys of `calib_funcs` as channel numbers. `datastore` must be a
(typically disk-based) dict of "ch..." channel names and raw data tables.
The data is read and filtered in chunks of size `chunk_size` internally.

Returns a dict by channel over dicts of raw data tables with the same
keys as `energy_windows`.
"""
function filter_raw_data_by_energy(
    datastore::AbstractDict{<:AbstractString,<:Any},
    calib_funcs::IdDict{<:Integer,<:Function},
    energy_windows::AbstractDict{Symbol,<:AbstractInterval{<:Number}};
    chunk_size = 10000
)
    channels = sort(collect(keys(calib_funcs)))
    IdDict((
        ch => map_chunked(get_raw_ch_data(datastore, ch), chunk_size) do chunk
            filter_raw_data_by_energy(chunk, calib_funcs[ch], energy_windows)
        end for ch in channels
    ));
end
=#
