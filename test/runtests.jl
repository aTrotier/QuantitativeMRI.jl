using qMRI
using Test

@testset "qMRI.jl" begin
    t = 1:1:100
    T2 = 10 # ms
    M0 = 1

    S = M0 * exp.(-t/T2)
    S = reshape(S,1,1,1,:)

    M0_fit,T2_fit = T2ExponentialFit(S,t)

    @test abs.(M0 - M0_fit[1]) < 10e-6 && abs.(T2 - T2_fit[1]) < 10e-6
end
