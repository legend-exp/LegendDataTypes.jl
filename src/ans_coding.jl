# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).


_binary_log(n::T) where {T<:Integer} = 8*sizeof(T) - leading_zeros(n+n-one(T)) - 1


function _ans_gen_pdf_cdf(::Type{TData}, symbol_pmf::Function) where TData<:Real
    TStream = widen(unsigned(TData)) # TBuf = widen(TStream) is implied
    m = Int(typemax(TStream)) - Int(typemin(TStream)) + 1
    n = _binary_log(m)
    @assert m == 1<<n

    n_symbols = Int(typemax(TData)) - Int(typemin(TData)) + 1

    s_cutoff::Int = 0
    total_prob::Float64 = 0
    while total_prob < 0.99865 && s_cutoff < (n_symbols - 1)
        total_prob += symbol_pmf(s_cutoff)
        s_cutoff += 1
    end

    float_PDF = [symbol_pmf(s) for s in 0:s_cutoff];
    # Tail entries have implicit PDF of one, add one to all explicit entries as well:
    float_PDF .= float_PDF .* ((m - n_symbols) / sum(float_PDF)) .+ 1

    PDF = round.(Int, float_PDF)
    n_tail_entries = n_symbols - length(eachindex(PDF))

    sum_PDF_err = (sum(PDF) + n_tail_entries) - m
    # ToDo: More efficient implementation:
    while sum_PDF_err != 0
        n_corr_idxs = min(abs(sum_PDF_err), length(eachindex(PDF)))
        err_sign = sign(sum_PDF_err)
        corr_idxs = sortperm(err_sign .* (float_PDF .- PDF))[begin:begin+n_corr_idxs-1]
        view(PDF, corr_idxs) .-= err_sign
        sum_PDF_err -= err_sign * length(eachindex(corr_idxs))
    end
    @assert sum(PDF) + n_tail_entries == m

    # In AND, CDF sums up PDF *until* current entry:
    CDF = circshift(cumsum(PDF), 1)
    CDF[begin] = 0

    (PDF = PDF, CDF = CDF, n_symbols = n_symbols)
end

function _get_ans_pdf(PDF::AbstractVector{<:Integer}, s::Integer, n_symbols::Integer)
    @assert 0 <= s < n_symbols
    T = eltype(PDF)
    @assert s < n_symbols
    s < length(eachindex(PDF)) ? PDF[firstindex(PDF) + s] : one(T)
end


function _get_ans_cdf(CDF::AbstractVector{<:Integer}, s::Integer, n_symbols::Integer)
    @assert 0 <= s < n_symbols
    T = eltype(CDF)
    l = length(eachindex(CDF))
    @assert s < n_symbols
    s < length(eachindex(CDF)) ? CDF[firstindex(CDF) + s] : T(last(CDF) + s - l + 1)
end


function _ans_encode!(
    output::AbstractVector{TStream}, data::AbstractVector{TData},
    PDF::AbstractVector{Int}, CDF::AbstractVector{Int}
) where {TStream <: Unsigned, TData <: Integer}
    TBuf = widen(TStream)

    m = Int(typemax(TStream)) - Int(typemin(TStream)) + 1
    n = _binary_log(m)
    n_symbols = Int(typemax(TData)) - Int(typemin(TData)) + 1

    d = (8 * sizeof(TBuf)) - n
    encstream_nbits = 8 * sizeof(TStream)

    @info "DEBUG A" TBuf m n n_symbols d encstream_nbits

    x::TBuf = zero(TBuf)
    for s in data
        pdf_s = TBuf(_get_ans_pdf(PDF, s, n_symbols))
        cdf_s = TBuf(_get_ans_cdf(CDF, s, n_symbols))
        @info "DEBUG B" x s pdf_s cdf_s

        if x >= (pdf_s << d)
            @info "DEBUG C push" unsafe_trunc(TStream, x)
            push!(output, unsafe_trunc(TStream, x))
            x = x >> encstream_nbits
        end

        # Accelerate by caching SignedMultiplicativeInverse{TBuf} for PDF:
        x = (div(x, pdf_s) << n) + mod(x, pdf_s) + cdf_s
    end

    while(x > 0)
        @info "DEBUG D push" x unsafe_trunc(TStream, x)
        push!(output, unsafe_trunc(TStream, x))
        x = x >> (8 * sizeof(TStream))
    end

    output
end


function _ans_decode!(data::AbstractVector{<:Integer}, input::AbstractVector{TStream}) where {TStream<:Unsigned}
    n = _binary_log(length(eachindex(PDF)))
    @assert axes(PDF) == axes(CDF)
    @assert sum(PDF) == 2^n #!!!!????

    TBuf = UInt64 #widen(TStream)
    mask = unsigned(2^n) - 1

    i0 = firstindex(PDF)
    i0 = firstindex(CDF)

    x::TBuf = zero(TBuf)

    nbits_input = 8 * sizeof(TStream)

    i_in = lastindex(input)
    for i_out in reverse(eachindex(data))
        while x < (1 << nbits_input) && i_in >= firstindex(input)
            x_low = input[i_in]
            i_in = i_in - 1
            x = (x << nbits_input) | x_low;
        end

        s = searchsortedlast(CDF, x & mask) - i0
        x = PDF[s + i0] * (x >> n) + (x & mask ) - CDF[s + i0]
    
        data[i_out] = s
    end

    data
end
