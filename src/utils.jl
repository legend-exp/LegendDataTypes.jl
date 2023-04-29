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

function _flatten_on_tables(xs)
    xs_1 = first(xs)
    ctor = Tables.materializer(xs_1)
    # ToDo: Avoid allcation due to broadcast of columns, if possible:
    ctor(_flatten_on_columns(Tables.columns.(xs)))
end

fast_flatten(xs::AbstractVector{<:StructArray{T}}) where T = StructArray{T}(_flatten_on_columns(StructArrays.components.(xs)))

# ToDo: Avoid copy due to map, if possible:
_append_lastdims_tplentry(tpls, ::Val{i}) where i = fast_flatten(map(x -> x[i], tpls))

@generated function _flatten_on_columns(nts::AbstractVector{<:NamedTuple{names}}) where names
    tpl_expr = :(())
    exprs = [:(_append_lastdims_tplentry(values(nts), Val($i))) for i in eachindex(names)]
    append!(tpl_expr.args, exprs)
    :(NamedTuple{$names}($tpl_expr))
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


"""
    map_chunked(f, table, chunk_size::Integer)

Maps a function `f` over a table in chunks of size `chunk_size`.

Calls `getindex` with contiguous index ranges, and so is also efficient for
disk-based arrays and similar arrays with slow serial indexing but fast
block-wise indexing.
"""
function map_chunked(f, table, chunk_size::Integer)
    idxs_partition = Iterators.partition(eachindex(table), chunk_size)
    fast_flatten([f(table[idxs]) for idxs in idxs_partition])
end
