export create_mp2rage

"""
    create_mp2rage(ima::Array{Complex{<:Real}})
    create_mp2rage(ima_magn::Array{<:Real},ima_phase::Array{<:Real})

Combine the image acquired at 2 differents inversion time (TI) in order to create a MP2RAGE image with the following equation.
ima could be of any size but the last dimension needs to be the 2 TI.

If magnitude and phase image are passed to the function. It will combine them to create a complex image.

``\\text{MP2RAGE} = Re(\\frac{ima_{TI_1} \\times conj(ima_{TI_2})}{|ima_{TI_1}|^2+|ima_{TI_2}|^2})``

# Arguments
- `ima::Array{Complex{<:Real}}`: image of any dimension with TI = 2
- `ima_magn::Array{<:Real},`: magnitude image of any dimension with TI = 2
- `ima_phase::Array{<:Real},`: magnitude image of any dimension with TI = 2

# Keywords

# Returns
- MP2RAGE images

# Bibliography
- Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele P-F, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage 2010;49:1271–1281 doi: 10.1016/j.neuroimage.2009.10.002.
"""
function create_mp2rage(ima::Array{Complex{T}}) where T<:Real
    ima_size = size(ima)
    @assert ima_size[end] == 2

    ima = reshape(ima,:,2)
    mp2 = real.(ima[:,1].*conj(ima[:,2])) ./ (abs.(ima[:,1]).^2 + abs.(ima[:,2]).^2 );
    mp2[isnan.(mp2)] .= 0
    mp2 = reshape(mp2,ima_size[1:end-1])
    return mp2
end

function create_mp2rage(ima_magn::Array{T},ima_phase::Array{T}) where T<:Real
    ima_size = size(ima_magn)
    phase_size = size(ima_phase)
    @assert ima_size == phase_size
    @assert ima_size[end]  == 2

    ima = ima_magn .* exp.(-im * ima_phase) # generate complex data
    return create_mp2rage(ima)
end
