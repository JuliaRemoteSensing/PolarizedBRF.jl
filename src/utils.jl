struct WignerDData{T}
    d₀₀::OffsetVector{T,Vector{T}}
    d₀₂::OffsetVector{T,Vector{T}}
    d₂₂::OffsetVector{T,Vector{T}}
    d₂₋₂::OffsetVector{T,Vector{T}}
    cjj₀₀::Vector{T}
    cjj₀₂::Vector{T}
    cjj₂₂::Vector{T}
    cj₀₀::Vector{T}
    cj₀₂::OffsetVector{T,Vector{T}}
    cj₂₂::OffsetVector{T,Vector{T}}
    j2::Vector{T}
    jj::Vector{T}
end

WignerDData(lmax::Integer) = begin
    lmax1 = lmax + 1
    d₀₀ = OffsetArray(zeros(lmax1), 0:lmax)
    d₀₂ = OffsetArray(zeros(lmax1), 0:lmax)
    d₂₂ = OffsetArray(zeros(lmax1), 0:lmax)
    d₂₋₂ = OffsetArray(zeros(lmax1), 0:lmax)
    cjj₀₀, cjj₀₂, cjj₂₂, cj₀₀, cj₀₂, cj₂₂, j2, jj = wignerd_coeff(lmax)
    return WignerDData(d₀₀, d₀₂, d₂₂, d₂₋₂, cjj₀₀, cjj₀₂, cjj₂₂, cj₀₀, cj₀₂, cj₂₂, j2, jj)
end

function wignerd_coeff(smax)
    cjj₀₀ = [Float64(j * (j + 1)^2) for j in 1:smax]
    cjj₀₂ = [j * (j + 1) * √((j + 1)^2 - 4) for j in 1:smax]
    cjj₂₂ = [Float64(j * ((j + 1)^2 - 4)) for j in 1:smax]
    cj₀₀ = [Float64((j + 1) * j^2) for j in 1:smax]
    cj₀₂ = OffsetArray([j * (j + 1) * √(j^2 - 4) for j in 2:smax], 2:smax)
    cj₂₂ = OffsetArray([Float64(j + 1) * (j^2 - 4) for j in 2:smax], 2:smax)
    j2 = [Float64(2j + 1) for j in 1:smax]
    jj = [Float64(j * (j + 1)) for j in 1:smax]

    return cjj₀₀, cjj₀₂, cjj₂₂, cj₀₀, cj₀₂, cj₂₂, j2, jj
end

function wignerd!(d, s, m, n, μ, cjj, cj, j2, jj)
    ss = max(abs(m), abs(n))
    c1 = ξ(m, n) * 2.0^(-ss) * √(factorial(2ss) / factorial(abs(m + n)) / factorial(abs(m - n)))

    d[ss] = c1 * (1 - μ)^(abs(m - n) / 2) * (1 + μ)^(abs(m + n) / 2)
    if ss == 0
        d[1] = μ
    end
    for j in max(1, ss):(s - 1)
        d[j + 1] = 1 / cjj[j] * (j2[j] * (jj[j] * μ - m * n) * d[j] - cj[j] * d[j - 1])
    end
end

"""
Calculate scattering matrix for given angles using the expansion coefficients.

$(SIGNATURES)

- `coeff`: Expansion coefficients, which is an `lmax1 x 6` matrix.
- `θ`: The angles in ascending order at which the scattering matrix is calculated.

Optional parameters:

- `output_dataframe`: Whether or not to output the scattering matrix in the `DataFrame` format. Default is `false`, and an `Nθ x 6` matrix will be outputed.
- `data`: Use precalculated Wigner d data. By default new data will be calculated.
"""
function F_from_expansion(coeff, θ; output_dataframe=false, data=WignerDData(size(coeff)[1] - 1))
    if θ[end] > π
        θ *= π / 180.0
    end

    Nθ = length(θ)
    μ = cos.(θ)
    F = zeros(Nθ, 6)
    lmax1 = size(coeff)[1]
    lmax = lmax1 - 1

    (; d₀₀, d₀₂, d₂₂, d₂₋₂, cjj₀₀, cjj₀₂, cjj₂₂, cj₀₀, cj₀₂, cj₂₂, j2, jj) = data

    for i in 1:Nθ
        μi = μ[i]
        wignerd!(d₀₀, lmax, 0, 0, μi, cjj₀₀, cj₀₀, j2, jj)
        wignerd!(d₀₂, lmax, 0, 2, μi, cjj₀₂, cj₀₂, j2, jj)
        wignerd!(d₂₂, lmax, 2, 2, μi, cjj₂₂, cj₂₂, j2, jj)
        wignerd!(d₂₋₂, lmax, 2, -2, μi, cjj₂₂, cj₂₂, j2, jj)

        F[i, 1] = sum(coeff[sj + 1, 1] * d₀₀[sj] for sj in 0:lmax)

        F₂₃ = sum((coeff[sj + 1, 2] + coeff[sj + 1, 3]) * d₂₂[sj] for sj in 0:lmax)
        F₂₋₃ = sum((coeff[sj + 1, 2] - coeff[sj + 1, 3]) * d₂₋₂[sj] for sj in 0:lmax)
        F[i, 2] = (F₂₃ + F₂₋₃) / 2
        F[i, 3] = (F₂₃ - F₂₋₃) / 2

        F[i, 4] = sum(coeff[sj + 1, 4] * d₀₀[sj] for sj in 0:lmax)

        # Note that Eqs. (7) & (8) in Mishchenko et al. (2016) do not have the 
        # negative sign, since associate Legendre function instead of Wigner d function
        # is used. But we need the negative sign here.
        F[i, 5] = -sum(coeff[sj + 1, 5] * d₀₂[sj] for sj in 0:lmax)
        F[i, 6] = -sum(coeff[sj + 1, 6] * d₀₂[sj] for sj in 0:lmax)
    end

    if output_dataframe
        df = DataFrame(; θ=θ * (180.0 / π), s11=F[:, 1], s22=F[:, 2], s33=F[:, 3], s44=F[:, 4], s12=F[:, 5],
                       s34=F[:, 6])
        return df
    else
        return F
    end
end
