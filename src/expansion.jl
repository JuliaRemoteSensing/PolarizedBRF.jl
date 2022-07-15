@enum ExpansionErrorType AbsoluteError RelativeError RootMeanSquaredError

@inline ξ(m, n) = n >= m ? 1 : (-1)^(m - n)

# TODO: find a better way to parrelize the expansion

function wignerd!(d, s, m, n, μ, cjj, cj, j2, jj)
    ss = max(abs(m), abs(n))
    c1 = ξ(m, n) * 2.0^(-ss) * √(factorial(2ss) / factorial(abs(m + n)) / factorial(abs(m - n)))

    d[ss] = c1 * (1 - μ)^(abs(m - n) / 2) * (1 + μ)^(abs(m + n) / 2)
    if ss == 0
        d[1] = μ
    end
    for j in max(1, ss):(s - 1)
        d[j + 1] = 1 / cjj[j] *
                   (j2[j] * (jj[j] * μ - m * n) * d[j] -
                    cj[j] * d[j - 1])
    end
end

"""
Expand given scattering matrix (6-independent columns: a1, a2, a3, a4, b1, b2) to General Spherical Function (GSF) coefficients.

$(SIGNATURES)

- `θ₀`: The scattering angles. If the maximum value is larger than `π`, then the angles are assumed to be in degrees, otherwise the angles are assumed to be in radians.
- `F`: The scattering matrix columns, which should be an `Nθ x 6` matrix.

Optional parameters:

- `smax`: Highest level for expansion. Default is `-1`, and in such case, the highest level is determined via trial-and-error until the desired precision (defined by `ε` and `error_type`) is reached.
- `smin`: The start point of iteration. Default is `30`.
- `double_threshold`: The highest level before which `s` is doubled each time. Default is `5000`.
- `adding_step`: The step size in the adding phase. Default is `100`.
- `stop_threshold`: The highest permitted expansion level. Default is `20000`.
- `ε`: The desired error bound. Default is `0.02`.
- `θₑ`: The angles used for error checking. Default is `θ₀`.
- `error_type`: The error type. Default is `PolarizedBRF.AbsoluteError`, other options include `PolarizedBRF.RelativeError` and `PolarizedBRF.RootMeanSquaredError`.
"""
function expand(θ₀, F;
                smax=-1, smin=30,
                double_threshold=5000,
                adding_step=100,
                stop_threshold=20000,
                ε=0.02, θₑ=collect(θ₀), error_type::ExpansionErrorType=AbsoluteError)
    if maximum(θ₀) > π
        θ₀ .*= π / 180.0
        θₑ .*= π / 180.0
    end

    itp = [LinearInterpolation(θ₀, F[:, i]; extrapolation_bc=Linear()) for i in 1:6]
    F11ₑ = itp[1].(θₑ)
    μₑ = cos.(θₑ)
    s = smax >= 0 ? smax : smin

    d₀₀ = OffsetArray(zeros(stop_threshold + 1), 0:stop_threshold)
    d₀₂ = OffsetArray(zeros(stop_threshold + 1), 0:stop_threshold)
    d₂₂ = OffsetArray(zeros(stop_threshold + 1), 0:stop_threshold)
    d₂₋₂ = OffsetArray(zeros(stop_threshold + 1), 0:stop_threshold)
    cjj₀₀ = [Float64(j * (j + 1)^2) for j in 1:stop_threshold]
    cjj₀₂ = [j * (j + 1) * √((j + 1)^2 - 4) for j in 1:stop_threshold]
    cjj₂₂ = [Float64(j * ((j + 1)^2 - 4)) for j in 1:stop_threshold]
    cj₀₀ = [Float64((j + 1) * j^2) for j in 1:stop_threshold]
    cj₀₂ = OffsetArray([j * (j + 1) * √(j^2 - 4) for j in 2:stop_threshold], 2:stop_threshold)
    cj₂₂ = OffsetArray([(j + 1) * (j^2 - 4) for j in 2:stop_threshold], 2:stop_threshold)
    j2 = [Float64(2j + 1) for j in 1:stop_threshold]
    jj = [Float64(j * (j + 1)) for j in 1:stop_threshold]

    while true
        μ, w = gausslegendre(s + 1)
        θ = acos.(μ)

        F_itp = [itp[i](θj) for i in 1:6, θj in θ]
        expansion = zeros(s + 1, 6)

        for i in 1:length(μ)
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

        F11ₑ′ = zeros(length(θₑ))
        for i in 1:length(θₑ)
            μi = μₑ[i]
            @views wignerd!(d₀₀, s, 0, 0, μi, cjj₀₀, cj₀₀, j2, jj)
            F11ₑ′[i] = sum(expansion[sj + 1, 1] * d₀₀[sj] for sj in 0:s)
        end

        maxerr = if smax < 0
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
            maxerr
        end

        if smax >= 0 || maxerr <= ε
            C = 1 / expansion[1, 1]
            expansion .*= C
            return expansion, C, maxerr
        end

        if s == stop_threshold
            @error "Expansion failed to reach desired precision."
            C = 1 / expansion[1, 1]
            expansion .*= C
            return expansion, C, maxerr
        end

        if s < double_threshold
            s = min(min(max(s, 1) * 2, double_threshold), stop_threshold)
        else
            s = min(s + adding_step, stop_threshold)
        end
    end
end
