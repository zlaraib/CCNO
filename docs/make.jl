using Documenter

# Automatically include all .jl files from the src directory
include("../src/CCNO.jl")

makedocs(
    sitename = "CCNO Documentation",
    format = Documenter.HTML(),
    modules = [CCNO],
    pages = [
        "Home" => "index.md"
    ],
    clean = true
)
