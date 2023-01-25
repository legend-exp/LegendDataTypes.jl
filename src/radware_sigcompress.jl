# This file is a part of LegendDataTypes.jl, licensed under the MIT License (MIT).

# The functions
#
# * _radware_sigcompress_encode_impl!
# * _radware_sigcompress_encode_impl!
#
# are almost literal translations of the functions compress_signal and
# decompress_signal from the radware-sigcompress v1.0 C-code (sigcompress.c
# as of 2017-10-25).
#
# The code was translated to Julia and is released here under MIT license with
# permission of the author of the original C-code
# (Copyright (c) 2018, David C. Radford <radforddc@ornl.gov>).


# C-style post-increment emulation
macro _pincr(x)
    quote
        begin
           tmp = $(esc(x))
           $(esc(x)) = $(esc(x)) + 1
           tmp
       end
    end
end


Base.@propagate_inbounds function _set_hton_u16!(A::AbstractVector{UInt8}, i, x::Integer)
    x_u16 = UInt16(x)
    i_1 = (i - firstindex(A)) * 2 + firstindex(A)
    i_2 = i_1 + 1
    A[i_1] = (x_u16 >> 8) % UInt8
    A[i_2] = (x_u16 >> 0) % UInt8
    return x
end


Base.@propagate_inbounds function _get_hton_u16(A::AbstractVector{UInt8}, i)
    i_1 = (i - firstindex(A)) * 2 + firstindex(A)
    i_2 = i_1 + 1
    UInt16(A[i_1]) << 8 | UInt16(A[i_2])
end


_get_high_u16(x::UInt32) = (x >> 16) % UInt16
_get_low_u16(x::UInt32) = (x >> 0) % UInt16

_set_high_u16(x::UInt32, y::UInt16) = x & UInt32(0x0000ffff) | (UInt32(y) << 16)
_set_low_u16(x::UInt32, y::UInt16) = x & UInt32(0xffff0000) | (UInt32(y) << 0)


const _radware_sigcompress_mask = UInt16[1, 3, 7,15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535]

function _radware_sigcompress_encode_impl!(sig_out::AbstractVector{UInt8}, sig_in::AbstractVector{<:Integer}, shift::Int32)
    @inbounds begin
        sig_len_in = length(eachindex(sig_in))

        max_output_len = 4 * sig_len_in
        if (length(eachindex(sig_out)) < max_output_len)
            resize!(sig_out, max_output_len)
        end

        i::Int32 = 0; j::Int32 = 0; max1::Int32 = 0; max2::Int32 = 0; min1::Int32 = 0; min2::Int32 = 0; ds::Int32 = 0; nb1::Int32 = 0; nb2::Int32 = 0
        iso::Int32 = 0; nw::Int32 = 0; bp::Int32 = 0; dd1::Int32 = 0; dd2::Int32 = 0

        dd::UInt32 = 0

        mask = _radware_sigcompress_mask

        #static int len[17] = {4096, 2048,512,256,128, 128,128,128,128,
        #                      128,128,128,128, 48,48,48,48};
        #= ------------ do compression of signal ------------ =#
        j = iso = 1; bp = 0

        _set_hton_u16!(sig_out, @_pincr(iso), UInt16(sig_len_in))  # signal length
        while (j <= sig_len_in)         # j = starting index of section of signal
            # find optimal method and length for compression of next section of signal
            si_j = Int16(sig_in[j] + shift)
            max1 = min1 = si_j
            max2 = Int32(-16000)
            min2 = Int32(16000)
            nb1 = nb2 = 2
            nw = 1
            i=j+1; while (i <= sig_len_in && i < j+48)  # FIXME; # 48 could be tuned better?
                si_i = Int16(sig_in[i] + shift)
                si_im1 = Int16(sig_in[i-1] + shift)
                if (max1 < si_i) max1 = si_i end
                if (min1 > si_i) min1 = si_i end
                ds = si_i - si_im1
                if (max2 < ds) max2 = ds end
                if (min2 > ds) min2 = ds end
                @_pincr(nw)
                @_pincr(i)
            end
            if (max1-min1 <= max2-min2)  # use absolute values
                nb2 = 99
                while (max1 - min1 > mask[nb1]) @_pincr(nb1) end
                #for (; i <= sig_len_in && i < j+len[nb1]; @_pincr(i)) {
                while (i <= sig_len_in && i < j+128)  # FIXME; # 128 could be tuned better?
                    si_i = Int16(sig_in[i] + shift)
                    if (max1 < si_i) max1 = si_i end
                    dd1 = max1 - min1
                    if (min1 > si_i) dd1 = max1 - si_i end
                    if (dd1 > mask[nb1]) break end
                    if (min1 > si_i) min1 = si_i end
                    @_pincr(nw)
                    @_pincr(i)
                end
            else                      # use difference values
                nb1 = 99
                while (max2 - min2 > mask[nb2]) @_pincr(nb2) end
                #for (; i <= sig_len_in && i < j+len[nb1]; @_pincr(i)) {
                while (i <= sig_len_in && i < j+128)  # FIXME; # 128 could be tuned better?
                    si_i = Int16(sig_in[i] + shift)
                    si_im1 = Int16(sig_in[i-1] + shift)
                    ds = si_i - si_im1
                    if (max2 < ds) max2 = ds end
                    dd2 = max2 - min2
                    if (min2 > ds) dd2 = max2 - ds end
                    if (dd2 > mask[nb2]) break end
                    if (min2 > ds) min2 = ds end
                    @_pincr(nw)
                    @_pincr(i)
                end
            end

            if (bp > 0) @_pincr(iso) end
            #=  -----  do actual compression  -----  =#
            _set_hton_u16!(sig_out, @_pincr(iso), nw)  # compressed signal data, first byte = # samples
            bp = 0               # bit pointer
            if (nb1 <= nb2)
                #=  -----  encode absolute values  -----  =#
                _set_hton_u16!(sig_out, @_pincr(iso), nb1)  # # bits used for encoding
                _set_hton_u16!(sig_out, @_pincr(iso), min1 % UInt16)  # min value used for encoding
                i = iso; while (i <= iso + div(nw*nb1, 16)) _set_hton_u16!(sig_out, i, 0); @_pincr(i) end
                i = j; while (i < j + nw)
                    si_i = Int16(sig_in[i] + shift)
                    dd = si_i - min1              # value to encode
                    dd = dd << (32 - bp - nb1)
                    _set_hton_u16!(sig_out, iso, _get_hton_u16(sig_out, iso) | _get_high_u16(dd))
                    bp += nb1
                    if (bp > 15)
                        _set_hton_u16!(sig_out, iso+=1, _get_low_u16(dd))
                        bp -= 16
                    end
                    @_pincr(i)
                end
            else
                #=  -----  encode derivative / difference values  -----  =#
                _set_hton_u16!(sig_out, @_pincr(iso), nb2 + 32)  # # bits used for encoding, plus flag
                _set_hton_u16!(sig_out, @_pincr(iso), si_j % UInt16)  # starting signal value
                _set_hton_u16!(sig_out, @_pincr(iso), min2 % UInt16)       # min value used for encoding
                i = iso; while (i <= iso + div(nw*nb2, 16)) _set_hton_u16!(sig_out, i, 0); @_pincr(i) end
                i = j+1; while (i < j + nw)
                    si_i = Int16(sig_in[i] + shift)
                    si_im1 = Int16(sig_in[i-1] + shift)
                    dd = si_i - si_im1 - min2     # value to encode
                    dd = dd << (32 - bp - nb2)
                    _set_hton_u16!(sig_out, iso, _get_hton_u16(sig_out, iso) | _get_high_u16(dd))
                    bp += nb2
                    if (bp > 15)
                        _set_hton_u16!(sig_out, iso+=1, _get_low_u16(dd))
                        bp -= 16
                    end
                    @_pincr(i)
                end
            end
            j += nw
        end

        if (bp > 0) @_pincr(iso) end
        if ((iso - 1)%2 > 0) @_pincr(iso) end     # make sure iso is even for 4-byte padding

        resize!(sig_out, 2 * (iso - 1))   # number of shorts in compressed signal data
    end

    return sig_out
end


function _radware_sigcompress_decode_impl!(sig_out::AbstractVector{<:Integer}, sig_in::AbstractVector{UInt8}, shift::Int32)
    @inbounds begin
        @assert iseven(length(eachindex(sig_in)))
        sig_len_in = div(length(eachindex(sig_in)), 2) # sig_in is UInt16 values in network byte order
        lastindex_sig_in = firstindex(sig_in) + div(lastindex(sig_in)-firstindex(sig_in), 2) # in terms of _get_hton_u16

        i::Int32 = 0; j::Int32 = 0; min::Int32 = 0; nb::Int32 = 0; isi::Int32 = 0; iso::Int32 = 0; nw::Int32 = 0; bp::Int32 = 0; siglen::Int32 = 0
        dd::UInt32 = 0
        mask = _radware_sigcompress_mask

        #= ------------ do decompression of signal ------------ =#
        j = isi = iso = 1; bp = 0
        siglen = _get_hton_u16(sig_in, @_pincr(isi)) % Int16  # signal length

        if (length(eachindex(sig_out)) < siglen)
            resize!(sig_out, siglen)
        end

        #printf("<<< siglen = %d\n", siglen);
        fill!(sig_out, zero(eltype(sig_out)))
        while (isi <= sig_len_in && iso <= siglen)
            if (bp > 0) @_pincr(isi) end
            bp = 0              # bit pointer
            nw = _get_hton_u16(sig_in,@_pincr(isi))  # number of samples encoded in this chunk
            nb = _get_hton_u16(sig_in,@_pincr(isi))  # number of bits used in compression

            if (nb < 32)
                #=  -----  decode absolute values  -----  =#
                min = _get_hton_u16(sig_in, @_pincr(isi)) % Int16  # min value used for encoding
                dd = _set_low_u16(dd, _get_hton_u16(sig_in, isi))
                i = 0; while (i < nw && iso <= siglen)
                    if (bp+nb > 15)
                        bp -= 16
                        dd = _set_high_u16(dd, _get_hton_u16(sig_in, @_pincr(isi)))
                        if isi <= lastindex_sig_in
                            dd = _set_low_u16(dd, _get_hton_u16(sig_in, isi))
                        end
                        dd = dd << (bp+nb)
                    else
                        dd = dd << nb
                    end
                    sig_out[@_pincr(iso)] = (_get_high_u16(dd) & mask[nb]) + min - shift
                    bp += nb
                    @_pincr(i)
                end
            else
                nb -= 32
                #=  -----  decode derivative / difference values  -----  =#
                sig_out[@_pincr(iso)] = _get_hton_u16(sig_in, @_pincr(isi)) % Int16 - shift  # starting signal value
                min = _get_hton_u16(sig_in, @_pincr(isi)) % Int16    # min value used for encoding
                if isi <= lastindex_sig_in
                    dd = _set_low_u16(dd, _get_hton_u16(sig_in, isi))
                end
                i = 1; while (i < nw && iso <= siglen)
                    if (bp+nb > 15)
                        bp -= 16
                        dd = _set_high_u16(dd, _get_hton_u16(sig_in, @_pincr(isi)))
                        if isi <= lastindex_sig_in
                            dd = _set_low_u16(dd, _get_hton_u16(sig_in, isi))
                        end
                        dd = dd << (bp+nb)
                    else
                        dd = dd << nb
                    end
                    sig_out[iso] = ((_get_high_u16(dd) & mask[nb]) + min + sig_out[iso-1] + shift) % Int16 - shift; @_pincr(iso)
                    bp += nb
                    @_pincr(i)
                end
            end
            j += nw
        end

        outlen = iso - 1
        if (siglen != outlen)
            throw(ErrorException("ERROR in decompress_signal: outlen $outlen != siglen $siglen"))
        end

        resize!(sig_out, siglen)
    end

    return sig_out
end



"""
    RadwareSigcompress <: AbstractArrayCodec
"""
struct RadwareSigcompress <: AbstractArrayCodec
    shift::Int
end
export RadwareSigcompress

RadwareSigcompress(::Type{T}) where {T<:Signed} = RadwareSigcompress(0)
RadwareSigcompress(::Type{T}) where {T<:Unsigned} = RadwareSigcompress(typemin(typeof(signed(zero(T)))))


function EncodedArrays.encode_data!(encoded::AbstractVector{UInt8}, codec::RadwareSigcompress, data::AbstractVector{T}) where {T}
    _radware_sigcompress_encode_impl!(encoded, data, Int32(codec.shift))
end


function EncodedArrays.decode_data!(data::AbstractVector{T}, codec::RadwareSigcompress, encoded::AbstractVector{UInt8}) where {T}
    _radware_sigcompress_decode_impl!(data, encoded, Int32(codec.shift))
end


function read_from_properties(read_property::Function, src::Any, ::Type{RadwareSigcompress})
    shft = read_property(src, :codec_shift, 0)
    RadwareSigcompress(shft)
end

function write_to_properties!(write_property::Function, dest::Any, codec::RadwareSigcompress)
    write_property(dest, :codec_shift, codec.shift)
    nothing
end
