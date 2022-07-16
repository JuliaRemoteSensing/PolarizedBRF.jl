using DelimitedFiles: readdlm
using PolarizedBRF
using Plots

coeff_raw = readdlm(joinpath(@__DIR__, "..", "fixture", "coeff2.txt"))
Cext = 2.84241
Csca = 2.57762
f = 0.2
reff = 0.6
λ = 0.63
lmax = size(coeff_raw)[1] - 1

(; coeff, F) = PolarizedBRF.ssf_correction(coeff_raw, Cext, Csca; volume_fraction=f, r=2π * reff / λ,
                                           interpolation_strategy=PolarizedBRF.Spline, strategy=PolarizedBRF.Ito18,
                                           output_dataframe=true)
plot(F.θ, F.s22 ./ F.s11; label="Ito18")

(; coeff, F) = PolarizedBRF.ssf_correction(coeff_raw, Cext, Csca; volume_fraction=f, r=2π * reff / λ,
                                           interpolation_strategy=PolarizedBRF.Spline, output_dataframe=true)
plot!(F.θ, F.s22 ./ F.s11; label="MAE <= 0.02")

(; coeff, F) = PolarizedBRF.ssf_correction(coeff_raw, Cext, Csca; volume_fraction=f, r=2π * reff / λ,
                                           interpolation_strategy=PolarizedBRF.Spline, Nm=4(size(coeff_raw)[1] - 1),
                                           output_dataframe=true)
plot(F.θ, F.s22 ./ F.s11; label="4lmax intermediate + MAE <= 0.02")

(; coeff, F) = PolarizedBRF.ssf_correction(coeff_raw, Cext, Csca; volume_fraction=f, r=2π * reff / λ,
                                           interpolation_strategy=PolarizedBRF.Spline, Nm=8(size(coeff_raw)[1] - 1),
                                           output_dataframe=true)
plot!(F.θ, F.s22 ./ F.s11; label="8lmax intermediate + MAE <= 0.02")

title!("F22/F11")
