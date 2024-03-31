using OAFCVX
using Documenter

DocMeta.setdocmeta!(OAFCVX, :DocTestSetup, :(using OAFCVX); recursive=true)

makedocs(;
    modules=[OAFCVX],
    authors="Mohammed Alshahrani <mmogib@gmail.com> and contributors",
    sitename="OAFCVX.jl",
    format=Documenter.HTML(;
        canonical="https://mmogib.github.io/OAFCVX.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mmogib/OAFCVX.jl",
    devbranch="master",
)
