using DelimitedFiles: readdlm
using PolarizedBRF

coeff = readdlm(joinpath(@__DIR__, "..", "fixture", "coeff.txt"))
NG = 49
x, R = PolarizedBRF.Wrapper.run_pbrf(1.0,
                                     NG,
                                     coeff;
                                     epsilon=1e-7,
                                     mode=PolarizedBRF.Standard)

@show R[:, :, 25, 25, 1]
@show R[:, :, 25, 25, 3]

itps = PolarizedBRF.get_itp(x, R)
RR = PolarizedBRF.evaluate(itps, 0.4, 0, 0.5, π)
@show 1.0 / π * 0.5 * RR * [1.0, 1.0, 0.0, 0.0]
