# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

function _ans_gen_pdf(::Type{TData}, symbol_pmf::Function) where TData<:Real
    TStream = widen(unsigned(TData)) # TBuf = widen(TStream) is implied
    m = typemax(TStream) - typemin(TStream) + 1
    n = _binary_log(typemax(TStream) - typemin(TStream) + 1)
    @assert m == 1<<n

    full_length_pdf = typemax(TData) - typemin(TData) + 1
    s_cutoff = min(ceil(Int, 2 * quantile(d_signed, cdf(Normal(), 3))), full_length_pdf - 1)

    float_PDF = [symbol_pmf(s) for s in 0:s_cutoff];
    # Tail entries have implicit PDF of one, add one to all explicit entries as well:
    float_PDF .= float_PDF .* ((m - full_length_pdf) / sum(float_PDF)) .+ 1

    PDF = round.(Int, float_PDF)
    n_tail_entries = full_length_pdf - length(eachindex(PDF))

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
end

function get_ans_pdf(PDF::AbstractVector{<:Integer}, s::Integer, m::Integer)
    @assert 0 <= s < m
    T = eltype(PDF)
    s < m ? PDF(firstindex(PDF) + s) : one(T)
end


function get_ans_cdf(CDF::AbstractVector{<:Integer}, s::Integer, m::Integer)
    @assert 0 <= s < m
    T = eltype(CDF)
    l = length(eachindex(CDF))
    s < m ? CDF(firstindex(CDF) + s) : T(last(CDF) + s - l + 1)
end


function _ans_encode!(output::AbstractVector{TStream}, data::AbstractVector{<:Integer}) where {TStream<:Unsigned}
    @assert axes(PDF) == axes(CDF)
    i0 = firstindex(PDF)

    TBuf = widen(TStream)
    sum_PDF = sum(PDF)
    n = _binary_log(sum_PDF)
    @assert sum_PDF == 1<<n

    d = (8 * sizeof(TBuf)) - n
    encstream_nbits = 8 * sizeof(TStream)

    x::TBuf = zero(TBuf)

    for s in data
        if x >= (PDF[s + i0] << d)
            push!(output, unsafe_trunc(TStream, x))
            x = x >> encstream_nbits
        end

        pdf_s = TBuf(PDF[s + i0])
        cdf_s = TBuf(CDF[s + i0])

        # Accelerate by caching SignedMultiplicativeInverse{TBuf} for PDF:
        x = (div(x, pdf_s) << n) + mod(x, pdf_s) + cdf_s
    end

    while(x > 0)
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
