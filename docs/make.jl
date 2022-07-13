using PolarizedBRF
using Documenter

DocMeta.setdocmeta!(PolarizedBRF, :DocTestSetup, :(using PolarizedBRF); recursive=true)

makedocs(;
    modules=[PolarizedBRF],
    authors="Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo="https://github.com/JuliaRemoteSensing/PolarizedBRF.jl/blob/{commit}{path}#{line}",
    sitename="PolarizedBRF.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaRemoteSensing.github.io/PolarizedBRF.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/JuliaRemoteSensing/PolarizedBRF.jl",
    devbranch="main"
)
