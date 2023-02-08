## ------------------------------------------------------------------
function _warmup(S, b, lb, ub, jump_args...)

    T = eltype(S)
    M, N = size(S)
    c = zeros(T, N)
    opm = FBAFluxOpModel(S, b, lb, ub, c, jump_args...)
    v0 = zeros(T, N)
    c0, c1 = zero(T), one(T)

    for i=1:N

        set_linear_obj!(opm, i, c1)
        optimize!(opm)
        v0 .+= solution(opm) / 2N
        
        set_linear_obj!(opm, i, -c1)
        optimize!(opm)
        v0 .+= solution(opm) / 2N

        set_linear_obj!(opm, i, c0)

    end
    return v0
end