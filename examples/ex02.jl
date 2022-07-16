using DelimitedFiles: readdlm
using PolarizedBRF

data = readdlm(joinpath(@__DIR__, "..", "fixture", "scattering_matrix.txt"))
θ₀ = data[:, 1]
F = data[:, 2:end]

R = PolarizedBRF.expand(F, θ₀; ε=0.02, error_type=PolarizedBRF.AbsoluteError)
