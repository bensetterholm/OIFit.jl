mas2rad(mas) = mas * pi / 648000000

function gauss(uv, inc, pa, fwhm, fluxg, fluxp)
    bl = baseline(uv, inc, pa)
    Θ = mas2rad(fwhm)
    return fluxg * exp.(-(π * Θ * bl).^2 / (4 * log(2))) .+ fluxp
end

function disk(uv, inc, pa, diam, fluxd, fluxp)
    bl = baseline(uv, inc, pa)
    Θ = mas2rad(diam)
    output = fill(1.0, size(bl))
    blmask = .!iszero.(bl)
    output[blmask] = 2 * besselj1.(π * Θ * bl[blmask]) ./ (π * Θ * bl[blmask])
    return fluxd * output .+ fluxp
end

function annulus(uv, inc, pa, diam, fluxa, fluxp)
    bl = baseline(uv, inc, pa)
    rho_0 = mas2rad(diam) / 2
    return fluxa * besselj0.(2 * π * rho_0 * bl) .+ fluxp
end

function pseudolorentzian(uv, inc, pa, hfdiam, fluxl, fluxp)
    bl = baseline(uv, inc, pa)
    theta = mas2rad(hfdiam)
    return fluxl * exp.(-π * theta * bl / sqrt(3)) .+ fluxp
end
