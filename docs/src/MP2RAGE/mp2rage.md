Here we will show how to reconstruct T1 maps from MP2RAGE images using the QuantitativeMRI.jl package.

## Step 1 : convert TI images to MP2RAGE / UNI image

We load the MP2RAGE dataset, which contains the two TI images.

```@example MP2RAGE
using QuantitativeMRI
using JLD2, CairoMakie
d = load(joinpath(pwd(),"../..","data","mp2rage.jld2"))

im_reco = d["im_reco"]

f=Figure(size=(800,600))
ax = Axis(f[1,1], title="TI1 image",aspect=1)
heatmap!(ax,abs.(im_reco[:,:,1]),colormap=:grays)
ax = Axis(f[1,2], title="TI2 image",aspect=1)
heatmap!(ax,abs.(im_reco[:,:,2]),colormap=:grays)
for ax in f.content
  hidedecorations!(ax)
end
f
```

We need to compute the MP2RAGE / UNI image with the function `mp2rage_comb`. This function expects to have the two TI images along the last dimension.

For example if you have multiple repetitions, you need to permute your images to have the TI images along the last dimension.

In our case, the images are in complex format but if you have the magnitude and phase images you can use the same function with `MP2 = mp2rage_comb(magnitude, phase)`.

```@example MP2RAGE
MP2 = mp2rage_comb(im_reco)
heatmap(MP2,colormap=:grays)
```

The range of data is now between -0.5 and 0.5.

## Step 2 : compute T1 maps
We can now compute the T1 maps using the function `mp2rage_T1maps`. This function expects the MP2RAGE image and the parameters of the MP2RAGE sequence.

The MP2RAGE parameters are stored in the `ParamsMP2RAGE` structure.

```raw
help?> ParamsMP2RAGE

  mutable struct ParamsMP2RAGE

  Fields
  ≡≡≡≡≡≡

  TI₁        :: Float64
  TI₂        :: Float64
  TR         :: Float64
  MP2RAGE_TR :: Float64
  ETL        :: Int64
  α₁         :: Float64
  α₂         :: Float64
```

```@example MP2RAGE
p = ParamsMP2RAGE(
  800,
  2250,
  6.5,
  5000,
  128,
  7,
  7
)
```

You can compute the corresponding look-up table (LUT) for cartesian encoding `mp2rage_lookuptable_cartesian(p::ParamsMP2RAGE;T1Range=1:0.5:10000,effInv = 0.96)` or radial`mp2rage_lookuptable_radial`.

Or directly compute the T1 maps as well as the LUT with `mp2rage_T1maps(MP2,p)`.
```@example MP2RAGE
T1, range_T1, LUT = mp2rage_T1maps(MP2,p,T1Range=1:0.5:5000,effInv = 0.96)

f=Figure()
ax = Axis(f[1,1], title="T1 map")
heatmap!(ax,T1,colorrange = (800,2000))
ax= Axis(f[1,2], title="LUT",xlabel="T1 [ms]", ylabel="MP2RAGE signal [a.u.]")
lines!(ax,range_T1,LUT)
f
```

