export T2Fit_Exp, T2Fit_ExpNoise, T2Fit_EpgNoise, qmri_echoAmplitudes

using MRISimulation
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
- fit_vec : fit objects for each pixels
- d : model dictionnary
    - d[:model]
    - d[:t]
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
    fit_vec = Vector{LsqFit.LsqFitResult}(undef,size(ima,2))

    @inbounds Threads.@threads for i in 1:size(ima, 2)
        try
            fit_vec[i] = curve_fit(model, t, ima[:, i], p0)
            fit_param[i,:] = fit_vec[i].param
        catch
            fit_param[i,:] .= NaN
        end
    end

    d=Dict{Symbol,Any}()
    d[:model]=model
    d[:t]=t

    return reshape(fit_param,dims[1:N-1]...,:),reshape(fit_vec,dims[1:N-1]...), d
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
- d : model dictionnary
    - d[:model]
    - d[:t]
    - d[:L] : number of coils used as optional parameters

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
    fit_vec = Vector{LsqFit.LsqFitResult}(undef,size(ima,2))

    @inbounds Threads.@threads for i = 1:size(ima, 2)
        try
            fit_vec[i] = curve_fit(model, t, ima[:, i], p0)
            fit_param[i,:] = fit_vec[i].param
        catch
            fit_param[i,:] .= NaN
        end
    end

    d=Dict{Symbol,Any}()
    d[:model]=model
    d[:t]=t
    d[:L]=L

    return reshape(fit_param,dims[1:N-1]...,:),reshape(fit_vec,dims[1:N-1]...), d
end

"""
    T2Fit_EpgNoise(ima::Array{T,4},params::Dict{Symbol,Any},x0::Vector{Float64}) where T<:Real

Fit the relaxation parameters T2 of a Multi-Spin Multi-Echo sequence with constant ``\\Delta TE`` and refocusing pulse.
EPG is computed with MRIReco functions.
In order to compute the gaussian noise standard deviation you need to know the number of coils `l`
and apply the following equation :

``\\sigma_g = \\frac{\\sigma}{\\sqrt{2L}}``

# Arguments
- `ima::Array{T,4}`: image with dimension (x,y,z,t)
- `t::Union{Vector{<:Real},StepRange{<:Real,<:Real}}`: times vector in ms
- `params::Dict{Symbol,Any}`:
- params[:TE]
- params[:train_length]
- params[:T1]

# Keywords


# Returns maps with format (x,y,z)
- M₀ maps (no unit)
- T₂ maps (ms)
- Delta maps (%)
- Noise maps (no unit)


# Bibliography
## Noise
- Cárdenas-Blanco A, Tejos C, Irarrazaval P, Cameron I. Noise in magnitude magnetic resonance images. Concepts Magn Reson Part A [Internet]. 2008 Nov;32A(6):409–16. Available from: http://doi.wiley.com/10.1002/cmr.a.20124
- Feng Y, He T, Gatehouse PD, Li X, Harith Alam M, Pennell DJ, et al. Improved MRI R 2 * relaxometry of iron-loaded liver with noise correction. Magn Reson Med [Internet]. 2013 Dec;70(6):1765–74. Available from: http://doi.wiley.com/10.1002/mrm.24607
## EPG
- MRIReco.jl implementation
"""
function T2Fit_EpgNoise(ima::Array{T,4}, params::Dict{Symbol,Any}, x0::Vector{Float64}=[1, 70.0, 0.9, 0]) where {T<:Real}
    params[:T1] = get(params, :T1, 1000.0)
    params[:delta] = get(params, :delta, 0.9)

    @assert haskey(params, :TE) && haskey(params, :train_length)
    @assert size(ima, 4) == params[:train_length]

    M0_maps = zeros(Float64, size(ima)[1:3])
    T2_maps = copy(M0_maps)
    delta_maps = copy(T2_maps)
    noise_maps = copy(T2_maps)

    for i in CartesianIndices(size(ima)[1:3])
        x0[1] = maximum(ima[i, :])
        fit_test = optimize(x -> residual(x, vec(ima[i, :]), params), x0)

        M0_maps[i] = fit_test.minimizer[1]
        T2_maps[i] = fit_test.minimizer[2]
        delta_maps[i] = fit_test.minimizer[3]
        noise_maps[i] = fit_test.minimizer[4]
    end

    return M0_maps, T2_maps, delta_maps, noise_maps
end


"""
qmri_echoAmplitudes(params::Dict{Symbol,Any})

Calculate EPG amplitudes of echos for a standard Multi-echo Spin-echo sequence
    (same TE / refoc pulse along the echo train):
- params[:TE]
- params[:train_length]
- params[:T1]
- params[:delta]
- params[:T2]
"""
function qmri_echoAmplitudes(params::Dict{Symbol,Any})
    TE = params[:TE]
    train_length = params[:train_length]
    delta = params[:delta]

    Refoc_vec = π * ones(Float64, train_length)
    TRefoc_vec = [TE/2:TE:TE*50...]
    TE_vec = [TE:TE:TE*train_length...]

    MESE_seq = MESequence(delta * π / 2, delta * Refoc_vec, TRefoc_vec, TE_vec)
    echos = echoAmplitudes(MESE_seq, 1 / params[:T1], 1 / params[:T2])
    return echos
end

"""
residual(x::Vector{Float64},ydata::Vector{Float64},params::Dict{Symbol,Any})

x -> Vector of parameter to fit :
    - x[1] : M0
    - x[2] : T2
    - x[3] : delta
    - x[4] : ``\\sigma``
"""
function residual(x::Vector{Float64}, ydata::Vector{Float64}, params::Dict{Symbol,Any})
    params[:T2] = x[2]
    params[:delta] = x[3]
    echos = x[1] * qmri_echoAmplitudes(params)
    # add noise
    echos_fit = sqrt.(abs.(echos) .^ 2 .+ x[4]^2)

    return sum((ydata - echos_fit) .^ 2)
end
