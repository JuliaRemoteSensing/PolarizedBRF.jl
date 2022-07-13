module PolarizedBRF

@enum QuadratureMode None NormalOriented Standard NeedNormal

using DocStringExtensions: SIGNATURES
using FastGaussQuadrature: gausslegendre
using Dierckx: Spline2D
using StaticArrays: @SMatrix

const Δ = @SMatrix [1 0 0 0
                    0 1 0 0
                    0 0 1 0
                    0 0 0 1]

const Δ₃₄ = @SMatrix [1 0 0 0
                      0 1 0 0
                      0 0 -1 0
                      0 0 0 -1]

const Δp = Δ + Δ₃₄
const Δm = Δ - Δ₃₄

include("wrapper.jl")

"""
Build interpolations from the calculated Fourier coefficients.

$(SIGNATURES)
"""
function get_itp(x, R)
    lmax1 = size(R)[5]

    return [Spline2D(x, x, R[i, j, :, :, m]) for i in 1:4, j in 1:4, m in 1:lmax1]
end


"""
Evaluate the refletion matrix ``\\mathbf{R}(\\mu,\\phi;\\mu_0,\\phi_0)`` using built interpolations.

$(SIGNATURES)
"""
function evaluate(itps, μ, ϕ, μ₀, ϕ₀)
    Δϕ = ϕ - ϕ₀
    mmax = size(itps)[3] - 1

    Rμϕ = @SMatrix(zeros(4, 4))
    for m in 0:mmax
        coeff = m == 0 ? 0.25 : 0.5
        Rm = @SMatrix([itps[i, j, m + 1](μ, μ₀) for i in 1:4, j in 1:4])
        Rμϕ += ((Δp * Rm * Δp + Δm * Rm * Δm) * cos(m * Δϕ) +
                (Δm * Rm * Δp - Δp * Rm * Δm) * sin(m * Δϕ)) * coeff
    end

    return Rμϕ
end

end
