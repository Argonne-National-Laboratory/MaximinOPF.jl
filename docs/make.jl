using Documenter

makedocs(    
    format = Documenter.HTML(analytics = "UA-367975-10", mathengine = Documenter.MathJax()),
    sitename = "CAPS",
    authors = "Sam",
    pages = [
        "Home" => "index.md"        
    ]
)

deploydocs(
    repo = "github.com/yim0331/MLSTOCK.git",
)