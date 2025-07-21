using QuantitativeMRI
using Documenter

DocMeta.setdocmeta!(QuantitativeMRI, :DocTestSetup, :(using QuantitativeMRI); recursive=true)

makedocs(;
    modules=[QuantitativeMRI],
    authors="Aurélien Trotier <a.trotier@gmail.com>",
    repo="https://github.com/aTrotier/QuantitativeMRI.jl/blob/{commit}{path}#{line}",
    sitename="QuantitativeMRI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aTrotier.github.io/QuantitativeMRI.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "MP2RAGE" => [
            "Standard reconstruction" => "MP2RAGE/mp2rage.md",
            "Slab profile correction" => "MP2RAGE/mp2rage_slice_profile.md",
            ],
    ],
)

deploydocs(;
    repo="github.com/aTrotier/QuantitativeMRI.jl",
    devbranch="main",
)
