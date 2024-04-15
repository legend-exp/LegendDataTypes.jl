# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendDataTypes

# Doctest setup
DocMeta.setdocmeta!(
    LegendDataTypes,
    :DocTestSetup,
    :(using LegendDataTypes);
    recursive=true,
)

makedocs(
    sitename = "LegendDataTypes",
    modules = [LegendDataTypes],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://LegendDataTypes.github.io/LegendDataTypes.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/LegendDataTypes.jl.git",
    forcepush = true,
    push_preview = true,
)
