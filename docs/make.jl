using SpacecraftBuilder
using Documenter

DocMeta.setdocmeta!(
    SpacecraftBuilder,
    :DocTestSetup,
    :(using SpacecraftBuilder);
    recursive = true,
)

makedocs(;
    modules = [SpacecraftBuilder],
    authors = "F. Capolupo",
    sitename = "SpacecraftBuilder.jl",
    format = Documenter.HTML(;
        canonical = "https://FraCpl.github.io/SpacecraftBuilder.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/FraCpl/SpacecraftBuilder.jl", devbranch = "master")
