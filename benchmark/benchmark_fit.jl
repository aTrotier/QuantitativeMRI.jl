using MRIReco
using qMRI
using Plots
using BenchmarkTools
plotly()

b = BrukerFile(pwd()*"/benchmark/data/MESE_FULLY")
ima = recoData(b)

TE = 7
t=TE:TE:TE*size(ima,4)

@time M0_fit,T2_fit = qMRI.T2ExpFit(ima[:,:,Int(end/2):Int(end/2),:],t);

@time M0_fit,T2_fit,σ_fit = qMRI.T2NoiseExpFit(ima[:,:,Int(end/2):Int(end/2),:],t,L=4); #48 seconds

heatmap(T2_fit[:,:,1],clims = (0,300),aspect_ratio = 1,legend = :none , axis=nothing)

