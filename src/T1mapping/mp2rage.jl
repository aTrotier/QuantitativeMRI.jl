export mp2rage_comb, ParamsMP2RAGE, mp2rage_T1maps, mp2rage_lookuptable

mutable struct ParamsMP2RAGE
    TI₁::Float64 #ms
    TI₂::Float64 #ms
    TR::Float64 #ms
    MP2RAGE_TR::Float64 #ms
    ETL::Int
    α₁::Float64 #degree
    α₂::Float64 #degree
end

"""
    mp2rage_comb(ima::Array{Complex{<:Real}})
    mp2rage_comb(ima_magn::Array{<:Real},ima_phase::Array{<:Real})

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
function mp2rage_comb(ima::Array{Complex{T}}) where T<:Real
    ima_size = size(ima)
    @assert ima_size[end] == 2

    ima = reshape(ima,:,2)
    mp2 = real.(ima[:,1].*conj(ima[:,2])) ./ (abs.(ima[:,1]).^2 + abs.(ima[:,2]).^2 );
    mp2[isnan.(mp2)] .= 0
    mp2 = reshape(mp2,ima_size[1:end-1])
    return mp2
end

function mp2rage_comb(ima_magn::Array{T},ima_phase::Array{T}) where T<:Real
    ima_size = size(ima_magn)
    phase_size = size(ima_phase)
    @assert ima_size == phase_size
    @assert ima_size[end]  == 2

    ima = ima_magn .* exp.(-im * ima_phase) # generate complex data
    return mp2rage_comb(ima)
end



"""
mp2rage_T1maps(im_MP2::Array{T},p::ParamsMP2RAGE;T1Range=1:10000,effInv = 0.96) where T <: Real

Compute Lookup table from MP2RAGE parameters

# Arguments
- `im_MP2::Array{T}`
- `p::ParamsMP2RAGE`


# Keywords
- T1Range = 1:10000
- effInv = 0.96

# Returns
- T1map
- T1Range
- lookUpTable

# Bibliography
- Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele P-F, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage 2010;49:1271–1281 doi: 10.1016/j.neuroimage.2009.10.002.

"""
function mp2rage_T1maps(im_MP2::Array{T},p::ParamsMP2RAGE;T1Range=1:10000,effInv = 0.96) where T<:Real
    # Generate lookUpTable + cut min
    lookUpTable, = mp2rage_lookuptable(p; T1Range=T1Range, effInv=effInv)
    maxVal,maxIdx = findmax(lookUpTable)
    T1Range = T1Range[maxIdx:end]
    lookUpTable = lookUpTable[maxIdx:end]

    minVal,minIdx = findmin(lookUpTable)
    T1Range = T1Range[1:minIdx]
    lookUpTable = lookUpTable[1:minIdx]

    T1map = MP2_T1.(im_MP2,Ref(lookUpTable),Ref(T1Range))
    return T1map, T1Range, lookUpTable
end

function MP2_T1(p_MP2,lookUpTable,T1Range)
    idxFirst = findfirst(lookUpTable .<= p_MP2)

    if !isnothing(idxFirst)
        T1 = T1Range[idxFirst]
    else
        T1 = 0
    end
    return T1
end

"""
    mp2rage_lookuptable(p::ParamsMP2RAGE;T1Range=1:0.5:10000,effInv = 0.96)

    Compute lookup table according to the MP2RAGE parameters

# Arguments
- `p::ParamsMP2RAGE`: MP2RAGE parameters structure

# Keywords
- `T1Range` : T1 range computed
- `effInv` : Inversion efficiency of the pulse

# Returns
- lookUpTable (NaN .= 0)

# Bibliography
- Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele P-F, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage 2010;49:1271–1281 doi: 10.1016/j.neuroimage.2009.10.002.
"""
function mp2rage_lookuptable(p::ParamsMP2RAGE;T1Range=1:0.5:10000,effInv = 0.96)
    TI1 =p.TI₁
    TI2 =p.TI₂
    TR = p.TR
    MP2RAGE_TR = p.MP2RAGE_TR
    n = p.ETL
    α₁ = p.α₁
    α₂ = p.α₂

    # compute exponential decay
    E1=vec(exp.(-TR./T1Range))
    EA=vec(exp.(-(TI1-(n./2-1).*TR)./T1Range))
    EB=vec(exp.(-(TI2-TI1-n.*TR)./T1Range))
    EC=vec(exp.(-(MP2RAGE_TR.-(TI2+(n./2).*TR))./T1Range))

    # compute mzss =[(B).*(cosd(α₂).*E1).^n+A].*EC+(1-EC);
    B = ((1 .- EA).*(cosd(α₁).*E1).^n + (1 .- E1).*(1 .- (cosd(α₁).*E1).^n)./(1 .- cosd(α₁).*E1)).*EB+(1 .- EB)
    A = (1 .- E1).*((1 .- (cosd(α₂).*E1).^n)./(1 .- cosd(α₂).*E1))

    mzss_num=((B).*(cosd(α₂).*E1).^n+A).*EC+(1 .-EC)
    mzss_denom=(1 .+ effInv.*(cosd(α₁).*cosd(α₂)).^n .* exp.(-MP2RAGE_TR./T1Range))
    mzss=mzss_num./mzss_denom

    # compute GRE1= sind(α₁).*(A.*(cosd(α₁).*E1).^(n./2-1)+B)
    A=-effInv.*mzss.*EA+(1 .- EA)
    B=(1 .- E1).*(1 .- (cosd(α₁).*E1).^(n./2-1))./(1 .- cosd(α₁).*E1)
    GRE1=sind(α₁).*(A.*(cosd(α₁).*E1).^(n./2-1)+B)

    # compute GRE2= sind(α₂).*(A-B)
    A=(mzss-(1 .- EC))./(EC.*(cosd(α₂).*E1).^(n./2))
    B=(1 .- E1).*((cosd(α₂).*E1).^(-n./2) .- 1)./(1 .- cosd(α₂).*E1)

    GRE2=sind(α₂).*(A-B)

    lookUpTable=GRE2.*GRE1./(GRE1.*GRE1+GRE2.*GRE2)
    lookUpTable[isnan.(lookUpTable)].= 0
    return lookUpTable, T1Range
end
