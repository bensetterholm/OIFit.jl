struct Vis2
    λ::Measurement{Float64}
    data::Vector{Measurement{Float64}}
    u::Vector{Float64}
    v::Vector{Float64}

    function Vis2(λ, data, u, v)
        if length(u) ≠ length(v)
            error("u and v arrays must have the same length")
        end
        if length(data) ≠ length(u)
            error("data array must have the same length as u and v arrays")
        end
        return new(λ, data, u, v)
    end

    function Vis2(
        λ::T,
        λerr::T,
        data::Vector{T},
        err::Vector{T},
        u::Vector{T},
        v::Vector{T},
    ) where T<:Real
        return Vis2(λ .± (0.5 .* λerr), data .± err, u, v)
    end
end

function ==(a::Vis2, b::Vis2)
    return a.λ == b.λ && a.data == b.data && a.u == b.u && a.v == b.v
end

function length(meas::Vis2)
    return length(meas.data)
end

"""
    purge(meas::Vis2)

Remove NaN values from a Vis2 object
"""
function purge(meas::Vis2)
    nan_val = isnan.(Measurements.value.(meas.data))
    nan_err = isnan.(Measurements.uncertainty.(meas.data))
    mask = .!(nan_val .| nan_err)
    return Vis2(meas.λ, meas.data[mask], meas.u[mask], meas.v[mask])
end

"""
    merge(arr::Vector{Vis2})

This merges Vis2 data with very similar eff_wave values. Currently, the implementation is
idiotic; I will need to think on this more in the future.
"""
function merge(arr::Vector{Vis2})
    λ = map(meas -> meas.λ, arr)
    m2c(x::Measurement{Float64}) = Complex(x.val, x.err)
    c2m(x::Complex{Float64}) = x.re ± x.im
    function bin(x::Vector{Measurement{Float64}})
        lo = minimum(map(y -> y.val - y.err, x))
        hi = maximum(map(y -> y.val + y.err, x))
        er = 0.5*(hi - lo)
        return Complex(mean(x).val, er)
    end
    d = Dict{Complex{Float64}, Vector{Int64}}()
    for i in eachindex(λ)
        test = false
        for k in keys(d)
            if abs(stdscore(c2m(k), λ[i])) < 0.5 * √(2)
                tmparr = vcat(pop!(d, k), i)
                d[bin(λ[tmparr])] = tmparr
                test = true
                break
            end
        end
        if test == false
            d[m2c(λ[i])] = [i]
        end
    end
    kd = sort(collect(keys(d)), by=real, rev=true)
    ret = Vector{Vis2}(undef, length(kd))
    for (i, el) in enumerate(kd)
        data = vcat([arr[j].data for j in d[el]]...)
        u = vcat([arr[j].u for j in d[el]]...)
        v = vcat([arr[j].v for j in d[el]]...)
        ret[i] = Vis2(c2m(el), data, u, v)
    end
    return ret
end

function +(a::Vis2, b::Vis2)
    ret = merge([a,b])
    if length(ret) > 1
        error("Wavelengths are incompatible")
    else
        return ret[1]
    end
end
