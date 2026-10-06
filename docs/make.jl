using Documenter, PiMTs

makedocs(;
    modules = [PiMTs],
    sitename = "PiMTs.jl",
    authors = "Markus Pirke",
    format = Documenter.HTML(;
        assets = ["assets/custom.css"],
        sidebar_sitename = true,
        collapselevel = 2,
        warn_outdated = true,
    ),
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Examples" => Any[
            "examples/charge_spectrum.md",
            "examples/fitting.md",
        ],
        "API" => "api.md"
    ],
    repo = Documenter.Remotes.URL(
        "https://git.ecap.work/xe91xote/PiMTs.jl/blob/{commit}{path}#L{line}",
        "https://git.ecap.work/xe91xote/PiMTs.jl"
    ),
)

deploydocs(;
  repo = "git.ecap.work/xe91xote/PiMTs.jl",
  devbranch = "main",
  push_preview=true
)
