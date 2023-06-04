module Wrapper

using DocStringExtensions: SIGNATURES
using Libdl: dlext, dlopen, dlsym
using ..PolarizedBRF

export run_pbrf

const LIBPBRF = joinpath(@__DIR__, "..", "shared", "libpbrf." * dlext)

"""
Calculate the Fourier coefficients of the refletion matrix using PolarizedBRF.

$(SIGNATURES)

- `ω` is the single scattering albedo, which should be within `(0, 1]`.
- `ngauss` determines the number of points used in integration.
- `coeff` is the expansion coefficients, which should be an `R x 6` matrix.

Optional parameters:

- `ε` is the threshold for convergence check. Detault is `1e-7`.
- `mode` determines the way to do the quadrature. Default is `PolarizedBRF.StandardQuadrature` (`NQUADR=2` in the original code), which is a normal Gaussian-Legendre quadrature within `[0, 1]`. Other options are `PolarizedBRF.NormalOrientedQuadrature` (`NQUADR=1` in the original code), `PolarizedBRF.ExplicitNormalQuadrature` (`NQUADR=3` in the original code) and `PolarizedBRF.CustomQuadrature`. Note that:
    - `PolarizedBRF.StandardQuadrature` might not be suitable if the refletion in the normal direction is needed.
    - For `PolarizedBRF.CustomQuadrature`, you need to input suitable nodes `x` and weights `w` yourself.
- `x` is the custom quadrature nodes within `[0, 1]`. Only used when `mode` is `PolarizedBRF.CustomQuadrature`.
- `w` is the corresponding quadrature weights. Only used when `mode` is `PolarizedBRF.CustomQuadrature`.

Results:

- `R` is a `4 x 4 x NG x NG x LMAX1` array, storing all the Fourier coefficients.
- `x` is the quadrature nodes.
- `w` is the quadrature weights.

"""
function run_pbrf(ω, ngauss, coeff; ε=1e-7, x=zeros(ngauss), w=zeros(ngauss),
                  mode::PolarizedBRF.QuadratureMode=PolarizedBRF.StandardQuadrature)
    0.0 < ω <= 1.0 || error("albedo must be within the range (0, 1]")
    lmax1, c = size(coeff)
    c == 6 || error("invalid expansion coefficients")
    R = zeros(4, 4, ngauss, ngauss, lmax1)

    dlopen(LIBPBRF) do libpbrf
        ccall(dlsym(libpbrf, :__pbrf_MOD_run_pbrf), Cvoid,
              (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
               Ptr{Float64},
               Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), ngauss,
              lmax1,
              Int(mode), ε, ω, coeff[:, 1], coeff[:, 2], coeff[:, 3], coeff[:, 4], coeff[:, 5], coeff[:, 6], x,
              w, R)
        return
    end

    return R, x, w
end

"""
Calculate ``d_{m0}^s(x)``, ``d_{m,-2}^s(x)`` and ``d_{m,2}^s(x)`` for ``s\\in[0,l_{\\max}]``.

$(SIGNATURES)
"""
function dd(x, lmax, m)
    lmax1 = lmax + 1
    d1 = zeros(lmax1)
    d2 = zeros(lmax1)
    d3 = zeros(lmax1)

    dlopen(LIBPBRF) do libpbrf
        ccall(dlsym(libpbrf, :__pbrf_MOD_dd), Cvoid,
              (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), x, lmax1, m, d1,
              d2, d3)
        return
    end

    return d1, d2, d3
end

end
