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

    function Vis2(λ::Measurement{T}, data::Measurement{T}, u::T, v::T) where T<:Real
        return Vis2(λ, [data], [u], [v])
    end
end

function ==(a::Vis2, b::Vis2)
    return a.λ == b.λ && a.data == b.data && a.u == b.u && a.v == b.v
end

function getindex(a::Vis2, @nospecialize vals)
    return Vis2(a.λ, a.data[vals], a.u[vals], a.v[vals])
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

function load(file::T) where T<:AbstractString
    f = OIFITS.load(file)
    db = OIFITS.select(f, "OI_VIS2")

    if length(db) ≠ 1
        error("Currently, only a single Vis2 table per oifits file is supported")
    end

    λ = db[1].ins[:eff_wave] .± (0.5 .* db[1].ins[:eff_band])
    data = db[1][:vis2data] .± db[1][:vis2err]
    u = db[1][:ucoord]
    v = db[1][:vcoord]

    return [Vis2(λ[i], data[i,:], u, v) |> purge for i in 1:length(λ)]
end

function load(files::Vector{T}) where T<:AbstractString
    return merge(vcat(load.(files)...))
end

function loaddir(dir::T) where T<:AbstractString
    fits = dir .* filter(x -> occursin(".fits", x), readdir(dir))
    oifits = dir .* filter(x -> occursin(".oifits", x), readdir(dir))
    files = vcat(fits, oifits)
    return load(files)
end

function uv(meas::Vis2)
    return (1/meas.λ.val) .* hcat(meas.u, meas.v)'
end

"""
    baseline(uv, inc=0, pa=0)

Converts a 2×n array of (u, v) points, or a Vis2 object, in into a baseline separation. If
inclination and position angle parameters are passed, the function computes the effective
baseline for a source stretched and rotated by that amount on the sky.
"""
function baseline(uv_arr::Array{Float64,2}, inc=0.0, pa=0.0)
    upvp = [
        (cosd(pa)*cosd(inc)) (-sind(pa)*cosd(inc))
        (sind(pa))           (cosd(pa))           
    ] * uv_arr
    return hypot.(upvp[1,:], upvp[2,:])
end

baseline(meas::Vis2, inc=0.0, pa=0.0) = baseline(uv(meas), inc, pa)

function setminerror(vis2, percent=5.0)
    for i in eachindex(vis2.data)
        err = vis2.data[i].val * percent / 100.
        if err > vis2.data[i].err
            vis2.data[i] = vis2.data[i].val ± err
        end
    end
end
