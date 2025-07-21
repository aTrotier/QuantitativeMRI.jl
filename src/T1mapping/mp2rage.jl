export mp2rage_comb, ParamsMP2RAGE, mp2rage_T1maps, mp2rage_lookuptable,
mp2rage_lookuptable_radial, mp2rage_lookuptable_cartesian

using Statistics
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

Compute Lookup table from MP2RAGE parameters. 
keywords radial = false means that only the central echo is used to compute the signal (standard method for cartesian acquisition)
If radial = true, all the echoes are sum.

# Arguments
- `im_MP2::Array{T}`
- `p::ParamsMP2RAGE`


# Keywords
- T1Range = 1:10000
- effInv = 0.96
- radial = false 

# Returns
- T1map
- T1Range
- lookUpTable

# Bibliography
- Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele P-F, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage 2010;49:1271–1281 doi: 10.1016/j.neuroimage.2009.10.002.
"""
function mp2rage_T1maps(im_MP2::Array{T},p::ParamsMP2RAGE;T1Range=1:10000,effInv = 0.96,radial=false) where T<:Real
    # Generate lookUpTable + cut min
    if !radial
        lookUpTable,T1Range = mp2rage_lookuptable_cartesian(p; T1Range=T1Range, effInv=effInv)
    else
        lookUpTable,T1Range = mp2rage_lookuptable_radial(p; T1Range=T1Range, effInv=effInv)
    end
    maxVal,maxIdx = findmax(lookUpTable)
    T1Range = T1Range[maxIdx:end]
    lookUpTable = lookUpTable[maxIdx:end]

    minVal,minIdx = findmin(lookUpTable)
    T1Range = T1Range[1:minIdx]
    lookUpTable = lookUpTable[1:minIdx]

    # reverse lookUpTable and Range for optimization purpose + convert to type of MP2
    T1Range = T.(T1Range)
    lookUpTable = T.(lookUpTable)
    T1map = MP2_T1.(im_MP2,Ref(lookUpTable[end:-1:1]),Ref(T1Range[end:-1:1]))
    return T1map, T1Range, lookUpTable
end

function MP2_T1(p_MP2::T,lookUpTable::Vector{T},T1Range::Vector{T}) where T
    idxFirst = searchsortedfirst(lookUpTable, p_MP2,lt= <)

    if idxFirst <= length(T1Range)-1
        T1 = T1Range[idxFirst]+(T1Range[idxFirst+1]-T1Range[idxFirst])*(p_MP2-lookUpTable[idxFirst])
    else
        T1 = T(0)
    end
    return T1
end

######################################
######## Convert T1 maps to MP2RAGE image 
######################################

"""
T1maps_MP2RAGE(T1map::Array{T},p::ParamsMP2RAGE;T1Range=1:10000,effInv = 0.96,radial=false) where T<:Real


Generates the MP2RAGE / UNI images from the T1 maps. Compute Lookup table from MP2RAGE parameters. 

keywords radial = false means that only the central echo is used to compute the signal (standard method for cartesian acquisition)
If radial = true, all the echoes are sum.

# Arguments
- `T1map::Array{T}`
- `p::ParamsMP2RAGE`


# Keywords
- T1Range = 1:10000
- effInv = 0.96
- radial = false 

# Returns
- MP2
- T1Range
- lookUpTable
"""
function T1maps_mp2rage(T1map::Array{T},p::ParamsMP2RAGE;T1Range=1:10000,effInv = 0.96,radial=false) where T<:Real
     # Generate lookUpTable + cut min
    if !radial
        lookUpTable,T1Range,LUT_TI1,LUT_TI2 = mp2rage_lookuptable_cartesian(p; T1Range=T1Range, effInv=effInv)
    else
        lookUpTable,T1Range,LUT_TI1,LUT_TI2 = mp2rage_lookuptable_radial(p; T1Range=T1Range, effInv=effInv)
    end
    maxVal,maxIdx = findmax(lookUpTable)
    T1Range = T1Range[maxIdx:end]
    lookUpTable = lookUpTable[maxIdx:end]
    LUT_TI1 = LUT_TI1[maxIdx:end]
    LUT_TI2 = LUT_TI2[maxIdx:end]

    minVal,minIdx = findmin(lookUpTable)
    T1Range = T1Range[1:minIdx]
    lookUpTable = lookUpTable[1:minIdx]
    LUT_TI1 = LUT_TI1[1:minIdx]
    LUT_TI2= LUT_TI2[1:minIdx]

    MP2, T1Range, lookUpTable = T1maps_mp2rage(T1map,lookUpTable, T1Range)
    TI1 = T1maps_TI(T1map,LUT_TI1, T1Range)
    TI2 = T1maps_TI(T1map,LUT_TI2, T1Range)
    return MP2, TI1, TI2, T1Range, lookUpTable, LUT_TI1, LUT_TI2
end



"""
T1maps_mp2rage(T1map::Array{T},lookUpTable::AbstractVector, T1Range::AbstractVector) where T<:Real


Generates the MP2RAGE / UNI images from the T1 maps. Compute Lookup table from MP2RAGE parameters. 

keywords radial = false means that only the central echo is used to compute the signal (standard method for cartesian acquisition)
If radial = true, all the echoes are sum.

# Arguments
- `T1map::Array{T}`
- lookUpTable::AbstractVector
- T1Range::AbstractVector


# Returns
- MP2
- T1Range
- lookUpTable
"""
function T1maps_mp2rage(T1map::Array{T},lookUpTable::AbstractVector, T1Range::AbstractVector) where T<:Real

    # reverse lookUpTable and Range for optimization purpose + convert to type of MP2
    T1Range = T.(T1Range)
    lookUpTable = T.(lookUpTable)
    MP2 = T1_MP2.(T1map,Ref(lookUpTable),Ref(T1Range))
    return MP2, T1Range, lookUpTable
end

function T1_MP2(T1map::T,lookUpTable::AbstractVector,T1Range::AbstractVector) where T
    idxFirst = searchsortedfirst(T1Range, T1map,lt= <=)

    if idxFirst <= length(T1Range)-1
        MP2 = lookUpTable[idxFirst]+(lookUpTable[idxFirst+1]-lookUpTable[idxFirst])*(T1map-T1Range[idxFirst])
    else
        MP2 = T(-0.5)
    end
    return MP2
end


## Generate TI images from T1 map

"""
T1maps_mp2rage(T1map::Array{T},TI_LUT::AbstractVector, T1Range::AbstractVector) where T<:Real


Generates the MP2RAGE / UNI images from the T1 maps. Compute Lookup table from MP2RAGE parameters. 

keywords radial = false means that only the central echo is used to compute the signal (standard method for cartesian acquisition)
If radial = true, all the echoes are sum.

# Arguments
- `T1map::Array{T}`
- TI_LUT::AbstractVector
- T1Range::AbstractVector


# Returns
- TI images
"""
function T1maps_TI(T1map::Array{T},TI_LUT::AbstractVector, T1Range::AbstractVector) where T<:Real
    T1Range = T.(T1Range)
    TI_LUT = T.(TI_LUT)
    TI = T1_TI.(T1map,Ref(TI_LUT),Ref(T1Range))
    return TI
end

function T1_TI(T1map::T,TI_LUT::AbstractVector,T1Range::AbstractVector) where T
    idxFirst = searchsortedfirst(T1Range, T1map,lt= <=)

    if idxFirst <= length(T1Range)-1
        TI = TI_LUT[idxFirst]+(TI_LUT[idxFirst+1]-TI_LUT[idxFirst])*(T1map-T1Range[idxFirst])
    else
        TI= T(0)
    end
    return TI
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
function mp2rage_lookuptable_cartesian(p::ParamsMP2RAGE;T1Range=1:0.5:10000,effInv = 0.96)
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
    return lookUpTable, T1Range, GRE1, GRE2
end

"""
    mp2rage_lookuptable_radial(p::ParamsMP2RAGE;T1Range=1:0.5:10000,effInv = 0.96)

    Compute lookup table according to the MP2RAGE parameters summing all the echoesl (mandatory for radial acquisition)

# Arguments
- `p::ParamsMP2RAGE`: MP2RAGE parameters structure

# Keywords
- `T1Range` : T1 range computed
- `effInv` : Inversion efficiency of the pulse

# Returns
- lookUpTable (NaN .= 0)

# Bibliography
- Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele P-F, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. NeuroImage 2010;49:1271–1281 doi: 10.1016/j.neuroimage.2009.10.002.
- Faller TL, Trotier AJ, Miraux S, Ribot EJ. Radial MP2RAGE sequence for rapid 3D T1 mapping of mouse abdomen: application to hepatic metastases. Eur Radiol. 2019 Nov;29(11):5844-5851. doi: 10.1007/s00330-019-06081-3. Epub 2019 Mar 19. PMID: 30888483.
"""
function mp2rage_lookuptable_radial(p::ParamsMP2RAGE;T1Range=1:0.5:10000,effInv = 0.96)
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

    GRE1_tmp=zeros(Float64,length(T1Range),n)
    GRE2_tmp=zeros(Float64,length(T1Range),n)
    for i in eachindex(1:n)
        # compute GRE1= sind(α₁).*(A.*(cosd(α₁).*E1).^(n./2-1)+B)
        A=-effInv.*mzss.*EA+(1 .- EA)
        B=(1 .- E1).*(1 .- (cosd(α₁).*E1).^(i-1))./(1 .- cosd(α₁).*E1)
        GRE1_tmp[:,i]=sind(α₁).*(A.*(cosd(α₁).*E1).^(i-1)+B)

        # compute GRE2= sind(α₂).*(A-B)
        A=(mzss-(1 .- EC))./(EC.*(cosd(α₂).*E1).^(n-i))
        B=(1 .- E1).*((cosd(α₂).*E1).^(-(n-i)) .- 1)./(1 .- cosd(α₂).*E1)

        GRE2_tmp[:,i]=sind(α₂).*(A-B)
    end

    GRE1 = mean(GRE1_tmp,dims=2)
    GRE2 = mean(GRE2_tmp,dims=2)

    lookUpTable=GRE2.*GRE1./(GRE1.*GRE1+GRE2.*GRE2)
    lookUpTable[isnan.(lookUpTable)].= 0
    return vec(lookUpTable), T1Range, vec(GRE1), vec(GRE2)
end