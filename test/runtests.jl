using DelimitedFiles: readdlm
using FastGaussQuadrature: gausslegendre
using PolarizedBRF
using Test

const BENCHMARK0 = [1.054548 -0.005217 0 0
                    -0.005217 0.262763 0 0
                    0 0 0.146687 0.032284
                    0 0 -0.032284 0.208790]

const BENCHMARK1 = [0.039998 0.007864 -0.014080 -0.000034
                    0.007864 0.045349 0.035117 -0.019115
                    0.014080 -0.035117 0.164078 -0.032249
                    -0.000034 -0.019115 0.032249 0.114666]

@testset "PolarizedBRF.jl" begin
    @testset "Benchmark results for PBRF calculation" begin
        coeff = readdlm(joinpath(@__DIR__, "..", "fixture", "coeff.txt"))
        NG = 49

        @testset "StandardQuadrature" begin
            R, x, w = PolarizedBRF.Wrapper.run_pbrf(1.0,
                                                    NG,
                                                    coeff;
                                                    ε=1e-7,
                                                    mode=PolarizedBRF.StandardQuadrature)

            @test isapprox(R[:, :, 25, 25, 1], BENCHMARK0;
                           atol=1e-6,
                           rtol=1e-6)

            @test isapprox(R[:, :, 25, 25, 2], BENCHMARK1;
                           atol=1e-5,
                           rtol=1e-5)
        end

        @testset "CustomQuadrature" begin
            x, w = gausslegendre(NG)
            @. x = (x + 1) * 0.5
            @. w = w * 0.5

            R, _ = PolarizedBRF.Wrapper.run_pbrf(1.0, NG, coeff;
                                                 ε=1e-7, mode=PolarizedBRF.CustomQuadrature, x=x, w=w)

            @test isapprox(R[:, :, 25, 25, 1], BENCHMARK0;
                           atol=1e-6,
                           rtol=1e-6)

            @test isapprox(R[:, :, 25, 25, 2], BENCHMARK1;
                           atol=1e-5,
                           rtol=1e-5)
        end
    end
end
