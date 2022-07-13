using DelimitedFiles: readdlm
using PolarizedBRF
using Test

@testset "PolarizedBRF.jl" begin
    @testset "Benchmark results" begin
        coeff = readdlm(joinpath(@__DIR__, "..", "fixture", "coeff.txt"))
        NG = 49
        x, R = PolarizedBRF.Wrapper.run_pbrf(1.0,
                                             NG,
                                             coeff;
                                             epsilon=1e-7,
                                             mode=PolarizedBRF.Standard)

        @test isapprox(R[:, :, 25, 25, 1],
                       [1.054548 -0.005217 0 0
                        -0.005217 0.262763 0 0
                        0 0 0.146687 0.032284
                        0 0 -0.032284 0.208790];
                       atol=1e-6,
                       rtol=1e-6)

        @test isapprox(R[:, :, 25, 25, 2],
                       [0.039998 0.007864 -0.014080 -0.000034
                        0.007864 0.045349 0.035117 -0.019115
                        0.014080 -0.035117 0.164078 -0.032249
                        -0.000034 -0.019115 0.032249 0.114666];
                       atol=1e-5,
                       rtol=1e-5)
    end
end
