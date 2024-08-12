group_by_evtno(table::TypedTables.Table) = TypedTables.Table(consgroupedview(table[sortperm(table.evtno)].evtno, Tables.columns(table)))


"""
    fast_flatten(vector_of_arrays)

Flattens a vector of arrays into a single array, by concatenating the
inner arrays along the last dimension.
"""
function fast_flatten end


fast_flatten(xs::AbstractVector{<:AbstractArray}) = _fast_flatten_generic(xs)

function _fast_flatten_generic(xs::AbstractVector{<:AbstractArray{T,N}}) where {T,N}
    x1 = first(xs)
    sz1 = size(x1)
    inner_sz = ntuple(i -> sz1[i], Val(N-1))
    lastdim_sz_sum = cumsum(size.(xs, N))
    pushfirst!(lastdim_sz_sum, 0)
    out_sz_N = lastdim_sz_sum[end]
    result = similar(x1, (inner_sz..., out_sz_N))
    firstidx_N = firstindex(result, N)
    idx_ranges = UnitRange.(firstidx_N .+ lastdim_sz_sum[begin:end-1], firstidx_N .+ lastdim_sz_sum[begin+1:end] .- 1)
    colons = ntuple(_ -> :, Val(N-1))
    # Can be parallelized:
    for i in 0:length(eachindex(xs))-1
        result[colons..., idx_ranges[firstindex(idx_ranges) + i]] = xs[firstindex(xs) + i]
    end
    return result
end

# ToDo: Improve implemenation:
fast_flatten(xs::AbstractVector{<:VectorOfArrays}) = vcat(xs...)

fast_flatten(xs::AbstractVector{<:VectorOfSimilarArrays{T,M}}) where {T,M} = VectorOfSimilarArrays{T,M}(fast_flatten(flatview.(xs)))

fast_flatten(xs::AbstractVector{<:AbstractArray{<:NamedTuple}}) = _flatten_on_tables(xs)

fast_flatten(xs::Base.ValueIterator) = fast_flatten(collect(xs))

function _flatten_on_tables(xs)
    xs_1 = first(xs)
    ctor = Tables.materializer(xs_1)
    # ToDo: Avoid allcation due to broadcast of columns, if possible:
    ctor(flatten_by_key(Tables.columns.(xs)))
end

fast_flatten(xs::AbstractVector{<:StructArray{T}}) where T = StructArray{T}(flatten_by_key(StructArrays.components.(xs)))



function fast_flatten(xs::AbstractVector{<:VectorOfEncodedArrays{T}}) where T
    @argcheck length(xs) >= 1
    codecs = (x -> x.codec).(xs)
    codec = first(codecs)
    @argcheck all(isequal(codec), codecs)
    innersizes = fast_flatten((x -> x.innersizes).(xs))
    encoded = fast_flatten((x -> x.encoded).(xs))
    return VectorOfEncodedArrays{T}(codec, innersizes, encoded)
end

function fast_flatten(xs::AbstractVector{<:VectorOfEncodedSimilarArrays{T}}) where T
    @argcheck length(xs) >= 1
    codecs = (x -> x.codec).(xs)
    codec = first(codecs)
    @argcheck all(isequal(codec), codecs)
    innersizes = (x -> x.innersize).(xs)
    innersize = first(innersizes)
    @argcheck all(isequal(innersize), innersizes)
    encoded = fast_flatten((x -> x.encoded).(xs))
    return VectorOfEncodedSimilarArrays{T}(codec, innersize, encoded)
end



"""
    flatten_by_key(data::AbstractVector{<:IdDict{<:Any, <:AbstractVector}})::IdDict

Flattens a vector of IdDicts into a single IdDict, by concatenating the
entries for each key separately.
"""
function flatten_by_key(data::AbstractVector{<:IdDict{<:Any, <:AbstractVector}})
    ks = keys(first(data))
    IdDict((k => fast_flatten([d[k] for d in data]) for k in ks))
end


# ToDo: Avoid copy due to map, if possible:
_append_lastdims_tplentry(tpls, ::Val{i}) where i = fast_flatten(map(x -> x[i], tpls))

@generated function flatten_by_key(nts::AbstractVector{<:NamedTuple{names}}) where names
    tpl_expr = :(())
    exprs = [:(_append_lastdims_tplentry(values(nts), Val($i))) for i in eachindex(names)]
    append!(tpl_expr.args, exprs)
    :(NamedTuple{$names}($tpl_expr))
end


_flatten_chunks(xs::AbstractVector{<:AbstractArray}) = fast_flatten(xs)
_flatten_chunks(xs::AbstractVector{<:IdDict}) = flatten_by_key(xs)

"""
    map_chunked(f, table, chunk_size::Integer)

Maps a function `f` over a table in chunks of size `chunk_size`.

Calls `getindex` with contiguous index ranges, and so is also efficient for
disk-based arrays and similar arrays with slow serial indexing but fast
block-wise indexing.
"""
function map_chunked(f, table, chunk_size::Integer)
    idxs_partition = Iterators.partition(eachindex(table), chunk_size)
    _flatten_chunks([f(table[idxs]) for idxs in idxs_partition])
end


"""
    decode_data(data)

Decode any encoded arrays present in (possibly nested) data.
"""
function decode_data end
export decode_data

decode_data(A::EncodedArray) = collect(A)
decode_data(A::VectorOfEncodedSimilarArrays) = collect.(A)
decode_data(A::AbstractVector{<:EncodedArray}) = collect.(A)

decode_data(data) = data
decode_data(x::Tuple) = map(decode_data, x)
decode_data(x::NamedTuple) = map(decode_data, x)
decode_data(tbl::AbstractVector{<:NamedTuple}) = similar_table(tbl, map(decode_data, Tables.columns(tbl)))
decode_data(wf::RDWaveform) = RDWaveform(decode_data(wf.time), decode_data(wf.signal))
decode_data(wfs::ArrayOfRDWaveforms) = ArrayOfRDWaveforms((decode_data(wfs.time), decode_data(wfs.signal)))


"""
    LegendDataTypes.similar_table(orig_tbl::AbstractVector{<:NamedTuple{names}}, cols::NamedTuple{names}) 

Return a table similar to `orig_tbl` from columns `cols`.
"""
function similar_table end

similar_table(orig_tbl::AbstractVector{<:NamedTuple{names}}, cols::NamedTuple{names}) where names = Tables.materializer(orig_tbl)(cols)

# Using Tables.materializer is not (always) type stable for StructArray, so specialize:
similar_table(::StructArray{T}, cols::NamedTuple{names}) where {names, T<:NamedTuple{names}} = StructArray{T}(cols)
