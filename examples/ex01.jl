using DelimitedFiles: readdlm
using PBRF

coeff = readdlm(joinpath(@__DIR__, "..", "fixture", "coeff.txt"))
NG = 49
x, R = PBRF.Wrapper.run_pbrf(1.0, NG, coeff; epsilon=1e-7, mode=PBRF.Standard)

@show R[:, :, 25, 25, 1]
@show R[:, :, 25, 25, 3]

itps = PBRF.get_itp(x, R)
RR = PBRF.evaluate(itps, 0.4, 0, 0.5, π)
@show 1.0 / π * 0.5 * RR * [1.0, 1.0, 0.0, 0.0]
