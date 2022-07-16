"""
Apply the static structure factor correction to the given expansion coeffcients.

$(SIGNATURES)

- `coeff`: The original expansion coeffcients, an `lmax1 x 6` matrix.
- `Cext`: Extinction cross section.
- `Csca`: Scattering cross section.

Optional parameters:

- `volume_fraction`: The volume fraction taken by the scatterers. Default is `0.2`.
- `r`: The volumetrically equivalent/effective radius of the scatterer.
- `Nθ`: Number of angles required for the output scattering matrix. Default is `181`.
- `Nm`: Number of angles used for the intermediate step. Default is `2lmax`.
- `equidistant`: Whether or not the output angles should be equidistant. Default is `true`. When set to `false`, Gaussian quadrature nodes (for cosθ) will be used instead.
- `output_dataframe`: Whether or not to output the scattering matrix in the `DataFrame` format. Default is `false`, and an `Nθ x 6` matrix will be outputed.
- `strategy`: Special strategy for expansion. Default is `PolarizedBRF.StandardExpansion`, and the expansion can be controlled via the keyword arguments specified for `expand()`. If use `PolarizedBRF.Ito18`, then the strategy used in Ito et al. (2018) will be used instead, in which the highest expansion level is set to `4 * lmax` and the expansion coefficients are cut off at the threshold `1e-8`.  
- All keyword arguments for `expand()` also apply here.
"""
function ssf_correction(coeff, Cext, Csca; volume_fraction=0.2, r=1.0, Nθ=181, Nm=2(size(coeff)[1] - 1),
                        equidistant=true, output_dataframe=false, strategy=StandardExpansion, kwargs...)
    lmax1 = size(coeff)[1]
    lmax = lmax1 - 1
    f = volume_fraction
    Cabs = Cext - Csca

    if equidistant
        θ = range(0, π, Nθ)
    else
        μ, _ = gausslegendre(Nθ)
        θ = reverse(acos.(μ))
    end

    dCsca = Csca / 4π
    Csca′ = 0.0

    μₘ, wₘ = gausslegendre(Nm)
    θₘ = reverse(acos.(μₘ))
    Fₘ = F_from_expansion(coeff, θₘ)

    for i in 1:length(θₘ)
        s = dCsca * ssf(f, r, θₘ[i])
        @views Fₘ[i, :] .*= s
        Csca′ += wₘ[i] * Fₘ[i, 1]
    end
    Csca′ *= 2π
    Cext′ = Csca′ + Cabs
    ω = Csca′ / Cext′
    P4CSca′ = 4π / Csca′
    Fₘ .*= P4CSca′

    # In Ito et al. (2018), number of quadrature points is set to `2lmax`,
    # and a cutoff threshold of `1e-8` is used, which is different from the 
    # practice of `spher_expan.f`.

    if strategy == Ito18
        coeff_new, _ = expand(Fₘ, θₘ; smax=4lmax, stop_threshold=4lmax, ngauss=_ -> 2lmax, cutoff=1e-8, kwargs...)
    else
        coeff_new, _ = expand(Fₘ, θₘ; kwargs...)
    end

    F = F_from_expansion(coeff_new, θ; output_dataframe=output_dataframe)

    return (coeff=coeff_new, F=F, Cext=Cext′, Csca=Csca′, Cabs=Cabs, g=coeff_new[2, 1] / 3.0, ω=ω)
end

"""
Calculate the static structure factor.

$(SIGNATURES)

- `f`: Volume fraction.
- `r`: Volumetrically equivalent radius.
- `θ`: Scattering angle in radians.
"""
function ssf(f, r, θ)
    u = 4r * sin(θ / 2)
    α = (1 + 2f)^2 / (1 - f)^4
    β = -6f * (1 + 0.5f)^2 / (1 - f)^4
    δ = 0.5α * f
    if abs(u) <= 0.01
        C = -24f * (α / 3 + β / 4 + δ / 6)
    else
        cosu = cos(u)
        sinu = sin(u)
        C = 24f * ((α + β + δ) / u^2 * cosu - (α + 2β + 4δ) / u^3 * sinu - (2β + 12δ) / u^4 * cosu +
             2β / u^4 +
             24δ / u^5 * sinu +
             24δ / u^6 * (cosu - 1))
    end

    return 1 / (1 - C)
end
