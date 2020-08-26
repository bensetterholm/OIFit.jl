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
