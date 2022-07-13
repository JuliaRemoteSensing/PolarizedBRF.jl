module Wrapper

using DocStringExtensions: SIGNATURES
using Libdl: dlopen, dlsym
using ..PBRF

const LIBPBRF = joinpath(@__DIR__, "..", "shared", "libpbrf" * (Sys.iswindows() ? ".dll" : Sys.islinux() ? ".so" : ".dylib"))

"""
Calculate the Fourier coefficients of the refletion matrix using PBRF.

$(SIGNATURES)

- `ω` is the single scattering albedo, which should be within `(0, 1]`.
- `ngauss` determines the number of points used in integration.
- `coeff` is the expansion coefficients, which should be an `R x 6` matrix.

Optional parameters:

- `epsilon` is the threshold for convergence check. Detault is `1e-7`.
- `mode` determines the way to do the quadrature. Default is `PBRF.Standard` (`NQUADR=2` in the original code), which is a normal Gaussian-Legendre quadrature within `[0, 1]`. Other options are `PBRF.NormalOriented` (`NQUADR=1` in the original code) and `PBRF.NeedNormal` (`NQUADR=3` in the original code). Note that `PBRF.Standard` might not be suitable if the refletion in the normal direction is needed.

Results:

- `x` is the quadrature nodes.
- `R` is a `4 x 4 x NG x NG x LMAX1` array, storing all the Fourier coefficients.
"""
function run_pbrf(ω, ngauss, coeff; epsilon=1e-7, mode::PBRF.QuadratureMode=PBRF.Standard)
    0.0 < ω <= 1.0 || error("albedo must be within the range (0, 1]")
    lmax1, c = size(coeff)
    c == 6 || error("invalid expansion coefficients")

    dlopen(LIBPBRF) do libpbrf
        x = zeros(ngauss)
        R = zeros(4, 4, ngauss, ngauss, lmax1)

        ccall(
            dlsym(libpbrf, :__pbrf_MOD_run_pbrf),
            Cvoid,
            (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
            ngauss,
            lmax1,
            Int(mode),
            epsilon,
            ω,
            coeff[:, 1],
            coeff[:, 2],
            coeff[:, 3],
            coeff[:, 4],
            coeff[:, 5],
            coeff[:, 6],
            x,
            R,
        )

        x, R
    end
end

end
