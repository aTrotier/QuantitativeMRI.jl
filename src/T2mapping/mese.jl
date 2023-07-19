export T2Fit_Exp, T2Fit_ExpNoise

"""
    T2Fit_Exp(ima::Array{T,N},t::AbstractVector{T},p0=nothing) where {T<:Real,N}

Fit the relaxation parameters T2 with the equation : ``S(t) = M_0 \\exp(-\\frac{t}{T2})``.

# Arguments
- `ima::Array{T,N}`: image with dimension [x,...,t]. Last dimensions -> temporal dimension
- `t::AbstractVector{<:Real}`: times vector in ms
- `p0=nothing`: starting values for fit, if empty p0=[maximum(ima),30]

# Keywords
- `removePoint::Bool=true`: remove the first point before fitting


# Returns
- fit_params : parameter maps
    - M₀ maps (no unit)
    - T₂ maps (ms)
"""
function T2Fit_Exp(ima::Array{T,N}, t::AbstractVector{T}, p0=nothing; removePoint::Bool=true) where {T<:Real,N}
    dims = size(ima)

    @assert dims[end] == length(t)
    ima = reshape(ima, :, dims[end])'

    if removePoint
        t = t[2:end]
        ima = ima[2:end,:]
    end

    if isnothing(p0)
        p0 = [maximum(ima), T.(30)]
    end

    model(t, p) = p[1] * exp.(-t / p[2])

    fit_param = zeros(eltype(ima),size(ima, 2),2)

    @inbounds Threads.@threads for i in 1:size(ima, 2)
        try
            fit_param[i,:] = curve_fit(model, t, ima[:, i], p0).param
        catch
            fit_param[i,:] .= NaN
        end
    end

    return reshape(fit_param,dims[1:N-1]...,:)
end

"""
    T2Fit_ExpNoise(ima::Array{T,N},t::AbstractVector{<:Real},p0=nothing; kwargs...) where {T<:Real,N}

Fit the relaxation parameters T2 with the equation : ``S(t) = \\sqrt{(M_0 \\exp(-\\frac{t}{T2}))^2 + 2 L \\sigma_g^2}``
where L est le nombre de canaux, et ``\\sigma_g`` le bruit gaussien sur les image

# Arguments
- `ima::Array{T,N}`: image with dimension [x,...,t]. Last dimensions -> temporal dimension
- `t::AbstractVector{<:Real}`: times vector in ms
- `p0=nothing`: starting values for fit, if empty p0=[maximum(ima),30,maximum(ima)*0.1]

# Keywords
- `removePoint::Bool=true`: remove the first point before fitting
- `L::Int=1`: Number of coil elements

# Returns
- fit_params : parameter maps
    - M₀ maps (no unit)
    - T₂ maps (ms)
    - Noise maps (no unit)
- fit_vec : fit objects for each pixels

# Bibliography
- Cárdenas-Blanco A, Tejos C, Irarrazaval P, Cameron I. Noise in magnitude magnetic
  resonance images. Concepts Magn Reson Part A [Internet]. 2008 Nov;32A(6):409–16. Available
  from: http://doi.wiley.com/10.1002/cmr.a.20124
- Feng Y, He T, Gatehouse PD, Li X, Harith Alam M, Pennell DJ, et al. Improved MRI R 2 *
  relaxometry of iron-loaded liver with noise correction. Magn Reson Med [Internet]. 2013
  Dec;70(6):1765–74. Available from: http://doi.wiley.com/10.1002/mrm.24607
"""
function T2Fit_ExpNoise(ima::Array{T,N}, t::AbstractVector{T}, p0=nothing; removePoint::Bool=true, L::Int=1) where {T<:Real,N}
    dims = size(ima)

    @assert dims[end] == length(t)
    ima = reshape(ima, :, dims[end])'

    if removePoint
        t = t[2:end]
        ima = ima[2:end,:]
    end

    if isnothing(p0)
        p0 = T.([maximum(ima), 30.0, maximum(ima) * 0.1])
    end

    model(t, p) = sqrt.((p[1] * exp.(-t / p[2])) .^ 2 .+ 2 * L * p[3]^2)

    fit_param = zeros(eltype(ima),size(ima, 2),3)
 
    @inbounds Threads.@threads for i = 1:size(ima, 2)
        try
            fit_param[i,:] = curve_fit(model, t, ima[:, i], p0).param
        catch
            fit_param[i,:] .= NaN
        end
    end

    return reshape(fit_param,dims[1:N-1]...,:)
end