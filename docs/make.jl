using Documenter
using JPEC

makedocs(
    sitename = "JPEC.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://OpenFUSIONToolkit.github.io/JPEC/"
    ),
    modules = [JPEC],
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "Splines" => "splines.md",
            "Vacuum" => "vacuum.md"
        ],
        "Examples" => [
            "Spline Examples" => "examples/splines.md",
            "Vacuum Examples" => "examples/vacuum.md",
            "Equilibrium Examples" => "examples/equilibrium.md"
        ]
    ],
    checkdocs = :exports
)

deploydocs(
    repo = "github.com/OpenFUSIONToolkit/JPEC.git",
    branch = "gh-pages",
    devbranch = "main"
)
