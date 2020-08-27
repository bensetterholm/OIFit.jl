const hp = 1  # huber parameter
const php = 1 # pseudo huber parameter
const sp = 3  # saturation parameter

function l1norm(resid)
    return @. sign(resid) * sqrt(abs(resid))
end

function l2norm(resid)
    return resid
end

function huber(resid)
    return [abs(r) < hp ? r/√(2) : sign(r) * sqrt(hp * (abs(r) - hp/2)) for r in resid]
end

function pseudo_huber(resid)
    return @. sign(resid) * php * sqrt(sqrt(1 + (resid/php)^2) - 1)
end

function saturate(resid)
    return [abs(r) < sp ? r : sp * sign(r) for r in resid]
end

function makesensible(param::Vector{Float64})
    param[1] = (param[1] + 90) % 180
    if param[1] < 0
        param[1] += 180
    end
    param[1] = abs(param[1] - 90) # set inclination
    param[2] %= 180
    if param[2] < 0
        param[2] += 180
    end # set position angle
    for i in 1:length(param)÷3
        param[3*i] = abs(param[3*i])
    end
end

function fit(
    vis2::Vis2,
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
)
    uv = extract_uv(vis2)
    guessParam = [10., 10., 3., 0.5, 0.5]
    parinfo = CMPFit.Parinfo(5)
    parinfo[4].limited = (1, 1)
    parinfo[4].limits = (0, 1) # disk model flux
    parinfo[5].limited = (1, 1)
    parinfo[5].limits = (0, 1) # point source flux
    function cmpfit_callback(param::Vector{Float64})
        resid = (Measurements.value.(vis2.data) - model(uv, param...).^2) ./
            Measurements.uncertainty.(vis2.data)
        return statistic(resid)
    end
    cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
end

function fit(
    vis2::Vector{Vis2},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
)
    lengths = length.(vis2)
    len = length(lengths)
    splits = zeros(Int64, len+1)
    for i in eachindex(lengths)
        splits[i+1] = splits[i] + lengths[i]
    end
    uv = extract_uv.(vis2)
    guessParam = [10., 10., fill(0.5, 3*len)...]
    parinfo = CMPFit.Parinfo(length(guessParam))
    for i in 0:len-1
        guessParam[3+(3*i)] = 3.0
        parinfo[4+(3*i)].limited = (1, 1)
        parinfo[4+(3*i)].limits = (0, 1)
        parinfo[5+(3*i)].limited = (1, 1)
        parinfo[5+(3*i)].limits = (0, 1)
    end
    function cmpfit_callback(param::Vector{Float64})
        resid = vcat([(
            Measurement.value.(vis2[i].data) -
            model(uv[i], param[1], param[2], param[3*i], param[3*i+1], param[3*i+2]).^2
        ) ./ Measurement.uncertainty.(vis2[i].data) for i in 1:len]...)
        return statistic(resid)
    end
    cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
end

function sensiblefit(vis2::Union{Vis2, Vector{Vis2}}, model::Function; kwargs...)
    a = fit(vis2, model, kwargs...)
    makesensible(a.param)
    return a
end

function sizefit(
    vis2::Vector{Vis2},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
)
    lengths = length.(vis2)
    len = length(lengths)
    splits = zeros(Int64, len+1)
    for i in eachindex(lengths)
        splits[i+1] = splits[i] + lengths[i]
    end
    uv = extract_uv.(vis2)
    guessParam = [10., 10., 3.0, fill(0.5, 2*len)...]
    parinfo = CMPFit.Parinfo(length(guessParam))
    for i in 0:len-1
        parinfo[4+(2*i)].limited = (1, 1)
        parinfo[4+(2*i)].limits = (0, 1)
        parinfo[5+(2*i)].limited = (1, 1)
        parinfo[5+(2*i)].limits = (0, 1)
    end
    function cmpfit_callback(param::Vector{Float64})
        resid = vcat([(
            Measurements.value.(vis2[i].data) -
            model(uv[i], param[1], param[2], param[3], param[2*i+2], param[2*i+3]).^2
        ) ./ Measurements.uncertainty.(vis2[i].data) for i in 1:len]...)
        return statistic(resid)
    end
    cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
end

function randfit(
    vis2::Vector{Vis2},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
    itr = 10,
)
    lengths = length.(vis2)
    len = length(lengths)
    splits = zeros(Int64, len+1)
    for i in eachindex(lengths)
        splits[i+1] = splits[i] + lengths[i]
    end
    uv = extract_uv.(vis2)
    function cmpfit_callback(param::Vector{Float64})
        resid = vcat([(
            Measurements.value.(vis2[i].data) -
            model(uv[i], param[1], param[2], param[3*i], param[3*i+1], param[3*i+2]).^2
        ) ./ Measurements.uncertainty.(vis2[i].data) for i in 1:len]...)
        return statistic(resid)
    end
    function makefit()
        guessParam = [90., 180., fill(1.0, 3*len)...] .* rand(3*len+2)
        parinfo = CMPFit.Parinfo(length(guessParam))
        for i in 0:len-1
            guessParam[3+(3*i)] *= 6.0
            parinfo[4+(3*i)].limited = (1, 1)
            parinfo[4+(3*i)].limits = (0, 1)
            parinfo[5+(3*i)].limited = (1, 1)
            parinfo[5+(3*i)].limits = (0, 1)
        end
        return cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
    end
    a = [makefit() for i in 1:itr]
    b = a[argmin([a[i].bestnorm for i in 1:itr])]
    makesensible(b.param)
    return b
end

function randsfit(
    vis2::Vector{Vis2},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
    itr = 10,
)
    lengths = length.(vis2)
    len = length(lengths)
    splits = zeros(Int64, len+1)
    for i in eachindex(lengths)
        splits[i+1] = splits[i] + lengths[i]
    end
    uv = extract_uv.(vis2)
    function cmpfit_callback(param::Vector{Float64})
        resid = vcat([(
            Measurements.value(vis2[i].data) -
            model(uv[i], param[1], param[2], param[3], param[2*i+2], param[2*i+3]).^2
        ) ./ Measurements.uncertainty(vis2[i].data) for i in 1:len]...)
        return statistic(resid)
    end
    function makefit()
        guessParam = [90., 180., 6.0, fill(1.0, 2*len)...] .* rand(2*len+3)
        parinfo = CMPFit.Parinfo(length(guessParam))
        for i in 0:len-1
            parinfo[4+(2*i)].limited = (1, 1)
            parinfo[4+(2*i)].limits = (0, 1)
            parinfo[5+(2*i)].limited = (1, 1)
            parinfo[5+(2*i)].limits = (0, 1)
        end
        return cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
    end
    a = [makefit() for i in 1:itr]
    b = a[argmin([a[i].bestnorm for i in 1:itr])]
    makesensible(b.param)
    return b
end

function bootfit(
    vis2::Vector{Vis2},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
    itr = 10,
    straps = 10,
)
    t0 = time()
    lengths = length.(vis2)
    len = length(lengths)
    splits = zeros(Int64, len+1)
    for i in eachindex(lengths)
        splits[i+1] = splits[i] + lengths[i]
    end

    sp = Random.Sampler(Random.MersenneTwister, 1:len)
    c = fill(NaN, 9, straps)
    for strap in 1:straps
        vis2m = vis2[rand(sp, len)]
        uv = extract_uv.(vis2m)
        function cmpfit_callback(param::Vector{Float64})
            resid = vcat([(vis2m[i].data - model(uv[i], param[1], param[2], param[3], param[2*i+2], param[2*i+3]).^2) ./ vis2m[i].err for i in 1:len]...)
            return statistic(resid)
        end
        function makefit()
            guessParam = [90., 180., 6.0, fill(1.0, 2*len)...] .* rand(2*len+3)
            parinfo = CMPFit.Parinfo(length(guessParam))
            for i in 0:len-1
                parinfo[4+(2*i)].limited = (1, 1)
                parinfo[4+(2*i)].limits = (0, 1)
                parinfo[5+(2*i)].limited = (1, 1)
                parinfo[5+(2*i)].limits = (0, 1)
            end
            return cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
        end
        a = [makefit() for i in 1:itr]
        b = a[argmin([a[i].bestnorm for i in 1:itr])]
        makesensible(b.param)
        c[:,strap] = b.param
    end
    ret = randsfit(
        vis2,
        model,
        statistic = statistic,
        set_flux = set_flux,
        parinfo = parinfo,
        config = config,
        itr = itr*10
    )
    ret.perror = [std(c[i,:]) for i in 1:size(c)[1]]
    ret.covar = collect(transpose(c))
    ret.elapsed = time() - t0
    return ret
end
