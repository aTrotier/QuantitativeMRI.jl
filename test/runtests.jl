using qMRI
using Test

@testset "qMRI.jl" begin
    @info "test T2ExpFit"
    t = 1:1:100
    T2 = 10 # ms
    M0 = 1

    S = M0 * exp.(-t/T2)
    S = reshape(S,1,1,1,:)
    M0_fit,T2_fit = T2ExpFit(S,[t...])
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
end
