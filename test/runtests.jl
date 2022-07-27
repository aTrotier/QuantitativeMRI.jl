using qMRI
using Test
using Noise
using Random

@testset "qMRI.jl" begin
    Random.seed!(42)
    @info "test T2ExpFit"
    t = 1:1:100
    T2 = 10 # ms
    M0 = 1

    S = M0 * exp.(-t/T2)
    S = reshape(S,1,1,1,:)
    M0_fit,T2_fit = T2ExpFit(S,[t...])
    @test abs.(M0 - M0_fit[1]) < 10e-6 && abs.(T2 - T2_fit[1]) < 10e-6

    @info "test T2ExpFit + remove"
    t = 1:1:100
    T2 = 10 # ms
    M0 = 1

    S = M0 * exp.(-t/T2)
    S = reshape(S,1,1,1,:)
    M0_fit,T2_fit = T2ExpFit(S,[t...],removePoint=true)
    @test abs.(M0 - M0_fit[1]) < 10e-6 && abs.(T2 - T2_fit[1]) < 10e-6

    @info "test T2NoiseExpFit"
    t = 1:1:100
    T2 = 30 # ms
    M0 = 1000
    σ = 50
    L = 3

    S = sqrt.((M0 * exp.(-t/T2)).^2 .+ 2*L*σ^2)
    S = reshape(S,1,1,1,:)
    M0_fit,T2_fit,σ_fit = T2NoiseExpFit(S,t;L=L)
    @test abs.(M0 - M0_fit[1]) < 10e-6 && abs.(T2 - T2_fit[1]) < 10e-6 && abs.(σ - σ_fit[1]) < 10e-6

    # epg fit
    @info "fit epg"
    # static parameter
    params = Dict{Symbol,Any}()
    params[:T1] = 1000.0
    params[:train_length] = 50
    params[:TE] = 7.0#ms
    σ=0.1
    # moving
    params[:T2] = 45.0
    params[:delta] = 0.8# B1map

    echos = qmri_echoAmplitudes(params)
    ## add noise + M0
    M0 = 10
    echos_noise = abs.(add_gauss(M0*echos, σ))
    echos_noise = reshape(echos_noise,1,1,1,:)
    # fit and check results
    x0 = [1,70.0,0.9,0]
    M0_maps, T2_maps, delta_maps, noise_maps = T2EpgNoiseFit(echos_noise,params,x0)

    @test (abs(M0_maps[1]) - M0)/M0 < 0.05
    @test (abs(T2_maps[1]) - params[:T2])/params[:T2] < 0.05
    @test (abs(delta_maps[1]) - params[:delta])/params[:delta] < 0.05
    @test (abs(noise_maps[1])/sqrt(2) - σ)/σ < 0.05
end


