module qMRI

export T2ExpFit, T2NoiseExpFit
using LsqFit

"""
    T2ExpFit(ima::Array{T,4},t::Union{Vector{<:Real},StepRange{<:Real,<:Real}},p0=nothing) where T

Fit the relaxation parameters T2 with the equation : ``S(t) = M_0 \\exp(-\\frac{t}{T2}) ``.

# Arguments
- `ima::Array{T,4}`: image with dimension (x,y,z,t)
- `t::Union{Vector{<:Real},StepRange{<:Real,<:Real}}`: times vector in ms
- `p0=nothing`: starting values for fit, if empty p0=[maximum(ima),30]

# Keywords
- `removePoint::Bool=true`: remove the first point before fitting

# Returns maps with format (x,y,z)
- M₀ maps (no unit)
- T₂ maps (ms)

"""
function T2ExpFit(ima::Array{T,4},t::Union{Vector{<:Real},StepRange{<:Real,<:Real}},p0=nothing;removePoint::Bool=true) where T
    @assert size(ima,4) == length(t)
    if removePoint
        t=t[2:end]
        ima=ima[:,:,:,2:end]
    end

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

"""
    T2NoiseExpFit(ima::Array{T,4},t::Union{Vector{<:Real},StepRange{<:Real,<:Real}},p0=nothing; kwargs...) where T

Fit the relaxation parameters T2 with the equation : ``S(t) = \\sqrt{(M_0 \\exp(-\\frac{t}{T2}))^2 + 2 L \\sigma_g^2}``
where L est le nombre de canaux, et ``\\sigma_g`` le bruit gaussien sur les image

# Arguments
- `ima::Array{T,4}`: image with dimension (x,y,z,t)
- `t::Union{Vector{<:Real},StepRange{<:Real,<:Real}}`: times vector in ms
- `p0=nothing`: starting values for fit, if empty p0=[maximum(ima),30,maximum(ima)*0.1]

# Keywords
- `removePoint::Bool=true`: remove the first point before fitting
- `L::Int=1`: Number of coil elements

# Returns maps with format (x,y,z)
- M₀ maps (no unit)
- T₂ maps (ms)
- Noise maps (no unit)

# Bibliography
- Cárdenas-Blanco A, Tejos C, Irarrazaval P, Cameron I. Noise in magnitude magnetic resonance images. Concepts Magn Reson Part A [Internet]. 2008 Nov;32A(6):409–16. Available from: http://doi.wiley.com/10.1002/cmr.a.20124
- Feng Y, He T, Gatehouse PD, Li X, Harith Alam M, Pennell DJ, et al. Improved MRI R 2 * relaxometry of iron-loaded liver with noise correction. Magn Reson Med [Internet]. 2013 Dec;70(6):1765–74. Available from: http://doi.wiley.com/10.1002/mrm.24607
"""
function T2NoiseExpFit(ima::Array{T,4},t::Union{Vector{<:Real},StepRange{<:Real,<:Real}},p0=nothing;removePoint::Bool=true,L::Int=1) where T
    @assert size(ima,4) == length(t)

    if removePoint
        t=t[2:end]
        ima=ima[:,:,:,2:end]
    end

    if isnothing(p0); p0=[maximum(ima),30,maximum(ima)*0.1]; end
    model(t, p) = sqrt.((p[1] * exp.(-t / p[2])).^2 .+ 2*L*p[3]^2)
    
    #fit_vec = LsqFitResult[]
    M0 = zeros(Float64,size(ima)[1:3])
    T2 = zeros(Float64,size(ima)[1:3])
    noise = zeros(Float64,size(ima)[1:3])

    ima = reshape(ima,:,size(ima,4))

    for i = 1:size(ima,1)
        fit=curve_fit(model, t, ima[i,:], p0)
        M0[i] = fit.param[1]
        T2[i] = fit.param[2]
        noise[i] = fit.param[3]
    end

    return M0,T2,noise
end

end
