
function started_log(xs)
    xsnz = xs[xs .> 0.0]
    q1, q3 = quantile!(log(xsnz), [0.25, 0.75])
    iqr = q3 - q1
    cutoff = exp(q1 - 2.0 * iqr)
    c = minimum(xsnz[xsnz .> cutoff])
    return log(xs .+ c)
end

function started_log_cor(xs, ys)
    return cor(started_log(xs), started_log(ys))
end

