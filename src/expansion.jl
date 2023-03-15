@enum ExpansionErrorType AbsoluteError RelativeError RootMeanSquaredError
@enum ExpansionStrategy StandardExpansion Ito18
@enum InterpolationStrategy Linear Spline

@inline ξ(m, n) = n >= m ? 1 : (-1)^(m - n)

# TODO: find a better way to parallelize the expansion

"""
Expand given scattering matrix (6-independent columns: a1, a2, a3, a4, b1, b2) to General Spherical Function (GSF) coefficients.

$(SIGNATURES)

- `F`: The scattering matrix columns, which should be an `Nθ x 6` matrix.
- `θ₀`: The scattering angles. If the maximum value is larger than `π`, then the angles are assumed to be in degrees, otherwise the angles are assumed to be in radians.


Optional parameters:

- `smax`: Highest level for expansion. Default is `-1`, and in such case, the highest level is determined via trial-and-error until the desired precision (defined by `ε` and `error_type`) is reached.
- `smin`: The start point of iteration. Default is `30`.
- `ngauss`: The function for determine the number of Gaussian quadrature points. It takes the current expansion level plus 1, `s + 1`, as the input. Default is `identity`.
- `double_threshold`: The highest level before which `s` is doubled each time. Default is `5000`.
- `adding_step`: The step size in the adding phase. Default is `100`.
- `stop_threshold`: The highest permitted expansion level. Default is `20000`.
- `interpolation_strategy`: The interpolation strategy. Default is `PolarizedBRF.Linear`, which uses `LinearInterpolation` from `Interpolations.jl`. Optional is `PolarizedBRF.Spline`, which uses `Spline1D` from `Dierckx.jl`.
- `cutoff`: Cut off the expansion coefficients when the absolute value of `α₁` at this level is smaller than `cutoff`. Default is `0.0`, meaning no cutoff will be applied. Note that Ito et al. (2018) used `cutoff = 1e-8`.
- `ε`: The desired error bound. Default is `0.02`.
- `θₑ`: The angles used for error checking. Default is `θ₀`.
- `error_type`: The error type. Default is `PolarizedBRF.AbsoluteError`, other options include `PolarizedBRF.RelativeError` and `PolarizedBRF.RootMeanSquaredError`.
"""
function expand(F, θ₀; smax=-1, smin=30, ngauss=identity, double_threshold=5000, adding_step=100, stop_threshold=20000,
                cutoff=zero(Float64), interpolation_strategy::InterpolationStrategy=Linear, ε=0.02,
                θₑ=vcat([θ₀[1]], [(θ₀[i] + θ₀[i + 1]) / 2 for i in 1:(length(θ₀) - 1)], [θ₀[end]]),
                error_type::ExpansionErrorType=AbsoluteError, _kwargs...)
    if maximum(θ₀) > π
        θ₀ .*= π / 180.0
        θₑ .*= π / 180.0
    end

    if interpolation_strategy == Linear
        itp = [LinearInterpolation(θ₀, F[:, i]; extrapolation_bc=Line()) for i in 1:6]
    else # Spline
        itp = [Spline1D(θ₀, F[:, i]) for i in 1:6]
    end

    F11ₑ = itp[1].(θₑ)
    μₑ = cos.(θₑ)
    s = smax >= 0 ? smax : smin

    (; d₀₀, d₀₂, d₂₂, d₂₋₂, cjj₀₀, cjj₀₂, cjj₂₂, cj₀₀, cj₀₂, cj₂₂, j2, jj) = WignerDData(stop_threshold)

    while true
        μ, w = gausslegendre(ngauss(s + 1))
        θ = acos.(μ)

        F_itp = [itp[i](θj) for i in 1:6, θj in θ]
        expansion = zeros(s + 1, 6)

        for i in eachindex(μ)
            μi = μ[i]

            wignerd!(d₀₀, s, 0, 0, μi, cjj₀₀, cj₀₀, j2, jj)
            wignerd!(d₀₂, s, 0, 2, μi, cjj₀₂, cj₀₂, j2, jj)
            wignerd!(d₂₂, s, 2, 2, μi, cjj₂₂, cj₂₂, j2, jj)
            wignerd!(d₂₋₂, s, 2, -2, μi, cjj₂₂, cj₂₂, j2, jj)

            wi = w[i]
            for sj in 0:s
                expansion[sj + 1, 1] += (sj + 0.5) * wi * F_itp[1, i] * d₀₀[sj]
                expansion[sj + 1, 2] += (sj + 0.5) * wi * (F_itp[2, i] + F_itp[3, i]) * d₂₂[sj]
                expansion[sj + 1, 3] += (sj + 0.5) * wi * (F_itp[2, i] - F_itp[3, i]) * d₂₋₂[sj]
                expansion[sj + 1, 4] += (sj + 0.5) * wi * F_itp[4, i] * d₀₀[sj]
                expansion[sj + 1, 5] -= (sj + 0.5) * wi * F_itp[5, i] * d₀₂[sj]
                expansion[sj + 1, 6] -= (sj + 0.5) * wi * F_itp[6, i] * d₀₂[sj]
            end
        end

        for sj in 0:s
            expansion[sj + 1, 2], expansion[sj + 1, 3] = (expansion[sj + 1, 2] + expansion[sj + 1, 3]) / 2,
                                                         (expansion[sj + 1, 2] - expansion[sj + 1, 3]) / 2
        end

        s_valid = s

        # Cut-off strategy which is abandoned in `spher_expan.f` but used in `PackedMedia.f`
        if !iszero(cutoff)
            for sj in 3:s
                if abs(expansion[sj + 1, 1]) < cutoff
                    s_valid = sj
                    break
                end
            end
        end

        F11ₑ′ = zeros(length(θₑ))
        for i in 1:length(θₑ)
            μi = μₑ[i]
            wignerd!(d₀₀, s, 0, 0, μi, cjj₀₀, cj₀₀, j2, jj)
            F11ₑ′[i] = sum(expansion[sj + 1, 1] * d₀₀[sj] for sj in 0:s_valid)
        end

        if error_type == AbsoluteError
            maxerr = maximum(abs.(F11ₑ′ - F11ₑ))
        elseif error_type == RelativeError
            maxerr = 200.0 * maximum(abs(f1 - f2) / abs(f1 + f2) for (f1, f2) in zip(F11ₑ′, F11ₑ))
        else
            # Note that the formula used here is different from the original Fortran version.
            sqsum = sum((f1 - f2)^2 for (f1, f2) in zip(F11ₑ′, F11ₑ))
            maxerr = sqrt(sqsum / length(F11′))
        end
        @info "Smax = $s, max error = $maxerr"

        if smax >= 0 || maxerr <= ε
            C = 1 / expansion[1, 1]
            expansion .*= C
            return (expansion=expansion[1:(s_valid + 1), :], C=C, maxerr=maxerr)
        end

        if s == stop_threshold
            @error "Expansion failed to reach desired precision."
            C = 1 / expansion[1, 1]
            expansion .*= C
            return (expansion=expansion[1:(s_valid + 1), :], C=C, maxerr=maxerr)
        end

        if s < double_threshold
            s = min(max(s, 1) * 2, double_threshold, stop_threshold)
        else
            s = min(s + adding_step, stop_threshold)
        end
    end
end
