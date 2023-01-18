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


function _radware_sigcompress_encode_impl!(sig_out::AbstractVector{UInt16}, sig_in::AbstractVector{Int16})
    sig_len_in = length(eachindex(sig_in))

    max_output_len = 2 * sig_len_in
    if (length(eachindex(sig_out)) < max_output_len)
        resize!(sig_out, max_output_len)
    end

    i::Int32 = 0; j::Int32 = 0; max1::Int32 = 0; max2::Int32 = 0; min1::Int32 = 0; min2::Int32 = 0; ds::Int32 = 0; nb1::Int32 = 0; nb2::Int32 = 0
    iso::Int32 = 0; nw::Int32 = 0; bp::Int32 = 0; dd1::Int32 = 0; dd2::Int32 = 0

    db = UInt16[0, 0]
    dd = reinterpret(UInt32, db)

    mask = UInt16[1, 3, 7,15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535]


    #static int len[17] = {4096, 2048,512,256,128, 128,128,128,128,
    #                      128,128,128,128, 48,48,48,48};
    #= ------------ do compression of signal ------------ =#
    j = iso = 1; bp = 0

    sig_out[@_pincr(iso)] = sig_len_in     # signal length
    while (j <= sig_len_in)         # j = starting index of section of signal
        # find optimal method and length for compression of next section of signal 
        max1 = min1 = sig_in[j]
        max2 = -16000
        min2 = 16000
        nb1 = nb2 = 2
        nw = 1
        i=j+1; while (i <= sig_len_in && i < j+48)  # FIXME; # 48 could be tuned better?
            if (max1 < sig_in[i]) max1 = sig_in[i] end
            if (min1 > sig_in[i]) min1 = sig_in[i] end
            ds = Int32(sig_in[i]) - Int32(sig_in[i-1])
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
                if (max1 < sig_in[i]) max1 = sig_in[i] end
                dd1 = max1 - min1
                if (min1 > sig_in[i]) dd1 = max1 - sig_in[i] end
                if (dd1 > mask[nb1]) break end
                if (min1 > sig_in[i]) min1 = sig_in[i] end
                @_pincr(nw)
                @_pincr(i)
            end
        else                      # use difference values
            nb1 = 99
            while (max2 - min2 > mask[nb2]) @_pincr(nb2) end
            #for (; i <= sig_len_in && i < j+len[nb1]; @_pincr(i)) {
            while (i <= sig_len_in && i < j+128)  # FIXME; # 128 could be tuned better?
                ds = Int32(sig_in[i]) - Int32(sig_in[i-1])
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
        sig_out[@_pincr(iso)] = nw  # compressed signal data, first byte = # samples
        bp = 0               # bit pointer
        if (nb1 <= nb2)
            #=  -----  encode absolute values  -----  =#
            sig_out[@_pincr(iso)] = nb1                    # # bits used for encoding
            sig_out[@_pincr(iso)] = min1 % UInt16  # min value used for encoding
            i = iso; while (i <= iso + div(nw*nb1, 16)) sig_out[i] = 0; @_pincr(i) end
            i = j; while (i < j + nw)
                dd[1] = sig_in[i] - min1              # value to encode
                dd[1] = dd[1] << (32 - bp - nb1)
                sig_out[iso] |= db[2]
                bp += nb1
                if (bp > 15)
                    sig_out[iso+=1] = db[1]
                    bp -= 16
                end
                @_pincr(i)
            end
        else
            #=  -----  encode derivative / difference values  -----  =#
            sig_out[@_pincr(iso)] = nb2 + 32  # # bits used for encoding, plus flag
            sig_out[@_pincr(iso)] = sig_in[j] % UInt16  # starting signal value
            sig_out[@_pincr(iso)] = min2 % UInt16       # min value used for encoding
            i = iso; while (i <= iso + div(nw*nb2, 16)) sig_out[i] = 0; @_pincr(i) end
            i = j+1; while (i < j + nw)
                dd[1] = Int32(sig_in[i]) - Int32(sig_in[i-1]) - min2     # value to encode
                dd[1] = dd[1] << (32 - bp - nb2)
                sig_out[iso] |= db[2]
                bp += nb2
                if (bp > 15)
                    sig_out[iso+=1] = db[1]
                    bp -= 16
                end
                @_pincr(i)
            end
        end
        j += nw
    end

    if (bp > 0) @_pincr(iso) end
    if ((iso - 1)%2 > 0) @_pincr(iso) end     # make sure iso is even for 4-byte padding

    resize!(sig_out, iso - 1)   # number of shorts in compressed signal data
    sig_out
end #= compress_signal =#


function _radware_sigcompress_decode_impl!(sig_out::AbstractVector{Int16}, sig_in::AbstractVector{UInt16})
    sig_len_in = length(eachindex(sig_in))

    i::Int32 = 0; j::Int32 = 0; min::Int32 = 0; nb::Int32 = 0; isi::Int32 = 0; iso::Int32 = 0; nw::Int32 = 0; bp::Int32 = 0; siglen::Int32 = 0
    db = UInt16[0, 0]
    dd = reinterpret(UInt32, db)
    mask = UInt16[1, 3, 7,15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535]

    #= ------------ do decompression of signal ------------ =#
    j = isi = iso = 1; bp = 0
    siglen = sig_in[@_pincr(isi)] % Int16  # signal length

    if (length(eachindex(sig_out)) < siglen)
        resize!(sig_out, siglen)
    end

    #printf("<<< siglen = %d\n", siglen);
    fill!(sig_out, 0)
    while (isi <= sig_len_in && iso <= siglen)
        if (bp > 0) @_pincr(isi) end
        bp = 0              # bit pointer
        nw = sig_in[@_pincr(isi)]  # number of samples encoded in this chunk
        nb = sig_in[@_pincr(isi)]  # number of bits used in compression

        if (nb < 32)
            #=  -----  decode absolute values  -----  =#
            min = sig_in[@_pincr(isi)] % Int16  # min value used for encoding
            db[1] = sig_in[isi]
            i = 0; while (i < nw && iso <= siglen)
                if (bp+nb > 15)
                    bp -= 16
                    db[2] = sig_in[@_pincr(isi)]
                    if isi <= lastindex(sig_in)
                        db[1] = sig_in[isi]
                    end
                    dd[1] = dd[1] << (bp+nb)
                else
                    dd[1] = dd[1] << nb
                end
                sig_out[@_pincr(iso)] = (db[2] & mask[nb]) + min
                bp += nb
                @_pincr(i)
            end
        else
            nb -= 32
            #=  -----  decode derivative / difference values  -----  =#
            sig_out[@_pincr(iso)] = sig_in[@_pincr(isi)] % Int16  # starting signal value
            min = sig_in[@_pincr(isi)] % Int16    # min value used for encoding
            if isi <= lastindex(sig_in)
                db[1] = sig_in[isi]
            end
            i = 1; while (i < nw && iso <= siglen)
                if (bp+nb > 15)
                    bp -= 16
                    db[2] = sig_in[@_pincr(isi)]
                    if isi <= lastindex(sig_in)
                        db[1] = sig_in[isi]
                    end
                    dd[1] = dd[1] << (bp+nb)
                else
                    dd[1] = dd[1] << nb
                end
                sig_out[iso] = ((db[2] & mask[nb]) + min + sig_out[iso-1]) % Int16; @_pincr(iso)
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
    sig_out
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
    shift = convert(T, codec.shift)
    data_conv = Int16.(data .+ shift)
    output = similar(data_conv, UInt16)
    _radware_sigcompress_encode_impl!(output, data_conv)
    resize!(encoded, 0)
    append!(encoded, reinterpret(UInt8, hton.(output)))
    encoded
end


function EncodedArrays.decode_data!(data::AbstractVector{T}, codec::RadwareSigcompress, encoded::AbstractVector{UInt8}) where {T}
    encoded_conv = ntoh.(reinterpret(UInt16, encoded))
    output = similar(encoded_conv, Int16)
    tmp_data = _radware_sigcompress_decode_impl!(output, encoded_conv)
    copyto!(data, tmp_data)
    data .-= convert(T, codec.shift)

    data
end


function read_from_properties(read_property::Function, src::Any, ::Type{RadwareSigcompress})
    shft = read_property(src, :radware_sigcompress_shift, 0)
    RadwareSigcompress(shft)
end

function write_to_properties!(write_property::Function, dest::Any, codec::RadwareSigcompress)
    write_property(dest, :radware_sigcompress_shift, codec.shift)
    nothing
end
