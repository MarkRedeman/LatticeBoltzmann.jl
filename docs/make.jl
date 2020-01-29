using Documenter, LatticeBoltzmann

makedocs(
    modules = [LatticeBoltzmann],
    format = :html,
    checkdocs = :exports,
    sitename = "LatticeBoltzmann.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/MarkRedeman/LatticeBoltzmann.jl.git",
)
