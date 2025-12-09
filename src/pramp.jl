"""
     pramp(
        f;
        p = [0, 1],
        hmin = (p[end] - p[begin]) / 1000,
        hmax = (p[end] - p[begin]) / 10,
        h = hmax,
        hgrow = 1.2,
        hdegrow = 0.5,
        verbose = false
    )

Run 'f(p)'  successively with parameter values from the range given by 'p'.
Stepsize is given by `h`. If `f` doesn't throw an error, `h` is increased by the factor of 
`hgrow` unless it exceeds `hmax`,  and the next step is performed with the new value.

If `f(p)` throws an error, solution is retried with a lower value of `h=h*hdegrow`.
In this case, if `p==p[begin]` or `h<hmin`, the error is rethrown.
"""
function pramp(
        f;
        p = [0, 1],
        hmin = (p[end] - p[begin]) / 1000,
        hmax = p[end] - p[begin],
        h = hmax,
        hgrow = 1.2,
        hdegrow = 0.5,
        verbose = false
    )
    p[begin] >= 0.0 || throw(ArgumentError("pramp: parameter range should be positive"))
    p[end] >= p[begin] || throw(ArgumentError("pramp: parameter range should be monotone"))
    pcurrent = p[begin]
    first = true
    h = min(h, hmax)
    while (pcurrent < p[end]) || first
        try
            if first
                ptrial = pcurrent
            else
                ptrial = min(pcurrent + h, p[end])
            end
            f(ptrial)
            if verbose
                if first
                    println("pramp - success: p=$(pcurrent)")
                else
                    println("pramp - success: p=$(pcurrent) + $(h)")
                end
            end
            pcurrent = ptrial
            if first
                first = false
            else
                h = min(h * hgrow, hmax)
            end
        catch e
            if verbose
                if first
                    println("pramp - error: p=$(pcurrent)")
                else
                    println("pramp - error: p=$(pcurrent) + $(h)")
                end
            end
            h = h * hdegrow
            if h < hmin || first
                rethrow(e)
            end
        end
    end
    return
end
