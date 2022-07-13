using Libdl: dlopen, dlclose

cd(joinpath("..", "fortran"))

libpbrf = if Sys.iswindows()
    "libpbrf.dll"
elseif Sys.islinux()
    "libpbrf.so"
else
    "libpbrf.dylib"
end

run(`gfortran -O3 -shared -fPIC -funroll-loops -o $libpbrf pbrf.f90`)
if !isdir(joinpath("..", "shared"))
    mkdir(joinpath("..", "shared"))
end
mv(libpbrf, joinpath("..", "shared", libpbrf), force=true)
rm("pbrf.mod", force=true)

if Sys.iswindows()
    try
        dlclose(dlopen(joinpath("..", "shared", libpbrf)))
    catch e
        error("64-bit gfortran is required to build this package on Windows.")
    end
end
