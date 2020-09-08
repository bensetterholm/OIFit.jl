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

function makesensible!(param::Vector{Float64})
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
    vis2::Vector{Vis2},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    guess = nothing,
    parinfo = nothing,
    config = nothing,
    ritr = 0,
)
    lengths = length.(vis2)
    len = length(lengths)
    splits = zeros(Int64, len+1)
    for i in eachindex(lengths)
        splits[i+1] = splits[i] + lengths[i]
    end
    _uv = uv.(vis2)
    guessParam = if isnothing(guess)
        [45., 131.3, fill(0.5, 3*len)...]
    else
        s = 2 + 3 * len
        length(guess) ≠ s && error("`guess` should have $s parameters")
        guess
    end
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
            Measurements.value.(vis2[i].data) -
            model(_uv[i], param[1], param[2], param[3*i], param[3*i+1], param[3*i+2]).^2
        ) ./ Measurements.uncertainty.(vis2[i].data) for i in 1:len]...)
        return statistic(resid)
    end
    b = if ritr > 0
        function makefit()
            guessParam = [90, 180, vcat(fill([6, 1, 1], len)...)...] .* rand(3*len+2)
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
        a = [makefit() for i in 1:ritr]
        a[argmin([a[i].bestnorm for i in 1:ritr])]
    else
        cmpfit(cmpfit_callback, guessParam, parinfo=parinfo, config=config)
    end
    makesensible!(b.param)
    return b
end

function fit(
    vis2::Vis2,
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    guess = nothing,
    parinfo = nothing,
    config = nothing,
    ritr = 0,
)
    fit(
        [vis2],
        model;
        statistic = statistic,
        set_flux = set_flux,
        guess = guess,
        parinfo = parinfo,
        config = config,
        ritr = ritr,
    )
end

function bootfit(
    vis2::Union{Vis2, Vector{Vis2}},
    model::Function;
    statistic::Function = l1norm,
    set_flux = 0,
    parinfo = nothing,
    config = nothing,
    ritr = 10,
    straps = 10,
)
    t0 = time()
    _vis2 = vis2 isa Vis2 ? [vis2] : vis2
    len = length(_vis2)
    sp = [Random.Sampler(Random.MersenneTwister, 1:length(v)) for v in _vis2]
    params = fill(NaN, 3*len+2, straps)
    for strap in 1:straps
        params[:, strap] .= fit(
            [_vis2[i][rand(sp[i], length(_vis2[i]))] for i in eachindex(_vis2)],
            model;
            statistic=statistic,
            set_flux=set_flux,
            parinfo=parinfo,
            config=config,
            ritr=ritr,
        ).param
        makesensible!(params[:, strap])
    end
    ret = fit(
        _vis2,
        model,
        statistic = statistic,
        set_flux = set_flux,
        parinfo = parinfo,
        config = config,
        ritr = ritr*10
    )
    ret.perror = [std(params[i,:]) for i in 1:size(params)[1]]
    ret.covar = collect(transpose(params))
    ret.elapsed = time() - t0
    return ret
end
