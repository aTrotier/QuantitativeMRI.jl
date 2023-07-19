export MESE_EPG, T2Fit_EpgNoise
"""
    T2Fit_EpgNoise(ima::Array{T,4}, T1,TE,ETL,x0::Vector{Float64}) where T<:Real

Fit the relaxation parameters T2 of a Multi-Spin Multi-Echo sequence with constant ``\\Delta TE`` and refocusing pulse.
EPG is computed with MRIReco functions.
In order to compute the gaussian noise standard deviation you need to know the number of coils `l`
and apply the following equation :

``\\sigma_g = \\frac{\\sigma}{\\sqrt{2L}}``

# Arguments
- `ima::Array{T,N}`: image where the echoes are in the last dimension
- `t::Union{Vector{<:Real},StepRange{<:Real,<:Real}}`: times vector in ms
- `T1`: T1 relaxation time
- `TE`: Echo Time
- `ETL`: Echo Train Length

# Keywords

# Returns fitted maps with format (x,y,z, M0,T2,delta (b1), Noise)

# Bibliography
## Noise
- Cárdenas-Blanco A, Tejos C, Irarrazaval P, Cameron I. Noise in magnitude magnetic resonance images. Concepts Magn Reson Part A [Internet]. 2008 Nov;32A(6):409–16. Available from: http://doi.wiley.com/10.1002/cmr.a.20124
- Feng Y, He T, Gatehouse PD, Li X, Harith Alam M, Pennell DJ, et al. Improved MRI R 2 * relaxometry of iron-loaded liver with noise correction. Magn Reson Med [Internet]. 2013 Dec;70(6):1765–74. Available from: http://doi.wiley.com/10.1002/mrm.24607
## EPG
- MRIReco.jl implementation
"""
function T2Fit_EpgNoise(ima::Array{T,N}, T1,TE,ETL, x0::Vector{Float64}=[1, 70.0, 0.9, 0]) where {T<:Real,N}
  dims = size(ima)
  @assert dims[end] == ETL
  
  fit_param = zeros(eltype(ima),dims[1:end-1]...,4)
  @inbounds Threads.@threads for i in CartesianIndices(dims[1:end-1])
    x0[1] = maximum(ima[i, :])
    fit_test = optimize(x -> residual_EpgNoise(x, vec(ima[i, :]), T1,TE,ETL), x0)

    fit_param[i,:] = fit_test.minimizer
  end

  return fit_param
end

"""
    residual(x::Vector{<:Real}, ydata::Vector{<:Real},T1,TE,ETL)

x -> Vector of parameter to fit :
    - x[1] : M0
    - x[2] : T2
    - x[3] : delta
    - x[4] : ``\\sigma``
"""
function residual_EpgNoise(x::Vector{<:Real}, ydata::Vector{<:Real},T1,TE,ETL)
  M0,T2,delta,noise = x

  echos = M0 * MESE_EPG(T2,T1,delta,TE,ETL)
  # add noises
  echos_fit = sqrt.(abs.(echos) .^ 2 .+ noise^2)

  return sum((ydata - echos_fit) .^ 2)
end

"""
    MESE_EPG(T2,T1,delta,TE,ETL)

Calculate EPG amplitudes of echos for a standard Multi-echo Spin-echo sequence
    (same TE / refoc pulse along the echo train):
- T2
- T1
- delta : delta B1 [0-1]
- TE
- ETL : Echo train length
"""
function MESE_EPG(T2,T1,delta,TE,ETL)
  T = eltype(complex(T2))
  E = EPGStates([T(0.0)],[T(0.0)],[T(1.0)])
  echo_vec = Vector{Complex{eltype(T2)}}()

  E = epgRotation(E,pi/2*delta, pi/2)
  # loop over refocusing-pulses
  for i = 1:ETL
    E = epgDephasing(E,1)
    E = epgRelaxation(E,TE,T1,T2)
    E = epgRotation(E,pi*delta,0.0)
    E = epgDephasing(E,1)
    push!(echo_vec,E.Fp[1])
  end

  return (echo_vec)
end