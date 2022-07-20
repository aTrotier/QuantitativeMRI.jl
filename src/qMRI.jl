module qMRI

export T2ExponentialFit
using LsqFit

"""
    T2ExponentialFit(ima::Array{T,4},t=Any[],p0=nothing) where T

Fit the relaxation parameters T2 with the equation : ``S(t) = M_0 \\exp(-\\frac{t}{T2}) ``.

# Arguments
- `ima::Array{T,4}`: image with dimension (x,y,z,t)
- `t=Real[]`: times vector in ms
- `p0=nothing`: starting values for fit, if empty p0=[maximum(ima),30]

# Keywords
- `removePoint::Bool=true`: remove the first point before fitting

# Returns maps with format (x,y,z)
- M₀ maps (no unit)
- T₂ maps (ms)

"""
function T2ExponentialFit(ima::Array{T,4},t=Union(Real[],StepRange),p0=nothing;removePoint::Bool=true) where T
    @assert size(ima,4) == length(t)

    if isnothing(p0); p0=[maximum(ima),30]; end


    model(t, p) = p[1] * exp.(-t / p[2])
    
    #fit_vec = LsqFitResult[]
    M0 = zeros(Float64,size(ima)[1:3])
    T2 = zeros(Float64,size(ima)[1:3])

    ima = reshape(ima,:,size(ima,4))

    for i = 1:size(ima,1)
        fit=curve_fit(model, t, ima[i,:], p0)
        M0[i] = fit.param[1]
        T2[i] = fit.param[2]
    end

    return M0,T2
end

end
