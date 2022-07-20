using qMRI
using Documenter

DocMeta.setdocmeta!(qMRI, :DocTestSetup, :(using qMRI); recursive=true)

makedocs(;
    modules=[qMRI],
    authors="Aurélien Trotier <a.trotier@gmail.com>",
    repo="https://github.com/aTrotier/qMRI.jl/blob/{commit}{path}#{line}",
    sitename="qMRI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aTrotier.github.io/qMRI.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aTrotier/qMRI.jl",
    devbranch="main",
)
