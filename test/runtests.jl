using QuantitativeMRI
using Test
using Noise
using Random

@testset "QuantitativeMRI.jl" begin
    @testset "T2Fit" begin
        Random.seed!(42)
        @info "T2Fit"
        t = 1.0:1.0:100.0
        T2 = 10 # ms
        M0 = 1

        S = M0 * exp.(-t/T2)
        S = reshape(S,1,1,1,:)
        fit_param = T2Fit_Exp(S,eltype(S).(t))
        @test abs.(M0 - fit_param[1]) < 10e-6 && abs.(T2 - fit_param[2]) < 10e-6

        @info "T2Fit + remove"
        t = 1:1:100
        T2 = 10 # ms
        M0 = 1

        S = M0 * exp.(-t/T2)
        S = reshape(S,1,1,1,:)
        fit_param = T2Fit_Exp(S,eltype(S).(t),removePoint=true)
        @test abs.(M0 - fit_param[1]) < 10e-6 && abs.(T2 - fit_param[2]) < 10e-6

        @info "T2Fit_ExpNoise"
        t = collect(1.0:1.0:100.0)
        T2 = 30 # ms
        M0 = 1000
        σ = 50
        L = 3

        S = sqrt.((M0 * exp.(-t/T2)).^2 .+ 2*L*σ^2)
        S = reshape(S,1,1,1,:)
        fit_param = T2Fit_ExpNoise(S,eltype(S).(t);L=L)
        @test abs.(M0 - fit_param[1]) < 10e-6 && abs.(T2 - fit_param[2]) < 10e-6 && abs.(σ - fit_param[3]) < 10e-6
    end    
    @testset "T2Fit EPG" begin
        @info "T2Fit_EpgNoise"
        # static parameter
        T1 = 1000.0
        ETL = 50
        TE = 7.0#ms
        σ=0.1

        # moving
        T2 = 45.0
        delta = 0.8 # B1map
        M0 = 10.0

        echos = M0 * MESE_EPG(T2,T1,delta,TE,ETL)
        
        ## add noise
        echos_noise = abs.(add_gauss(echos, σ))

        # fit and check results
        x0 = [1,70.0,0.9,0]
        fit_param = T2Fit_EpgNoise(echos_noise,T1,TE,ETL,x0)

        @test (abs(fit_param[1]) - M0)/M0 < 0.05
        @test (abs(fit_param[2]) - T2)/T2 < 0.05
        @test (abs(fit_param[3]) - delta)/delta < 0.05
        @test (abs(fit_param[4])/sqrt(2) - σ)/σ < 0.05
    end

    @testset "MP2RAGE" begin
        # generate data
        sArray = [10,10,5,2]
        ima = ones(Complex,sArray...)
        ima[:,:,:,2] .= 2;

        phase = rand(Float64,sArray[1:3]...)*2π
        phase2 = repeat(phase,1,1,1,2)

        ima2 = ima .* exp.(im*phase)
        imMP = mp2rage_comb(ima2)
        @test all((imMP .- 0.4) .< eps(Float64))

        imMP2 = mp2rage_comb(abs.(ima2),phase2)
        @test all((imMP2 .- 0.4) .< eps(Float64))

        # T1 maps
        TI₁ = 800
        TI₂ = 2200
        TR = 7
        MP2RAGE_TR = 6250
        ETL = 128
        α₁ = 7
        α₂ = 7

        p = ParamsMP2RAGE(TI₁,TI₂,TR,MP2RAGE_TR,ETL,α₁,α₂)

        # test Broadcasting
        T = Float32
        T1map, = mp2rage_T1maps(T.([[0,0 ];;[ 0,0]]),p)
        @test any(abs.(T1map .- [[1396.0,1396.0];;[1396.0,1396.0]]) .< [[1.0,1.0];;[1.0,1.0]])
    end
end
