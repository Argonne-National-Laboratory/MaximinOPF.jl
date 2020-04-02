using Documenter

makedocs(    
    format = Documenter.HTML(analytics = "UA-367975-10", mathengine = Documenter.MathJax()),
    sitename = "MaximinOPF",
    authors = "Sam",
    pages = [
        "Home" => "index.md",
        "Manual" => [
        	"Getting Started" => "gettingstarted.md",
        	"API Document" => "API.md",
        	"Examples" => "example.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/yim0331/yim0331.github.io.git",
)