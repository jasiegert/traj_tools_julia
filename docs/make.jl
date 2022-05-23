using Documenter, TrajTools

makedocs(
    sitename="TrajTools Documentation",
    # prettyURLs disabled to make local browsing easier -> should be enabled for online hosting
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Reading trajectories" => "trajecio.md",
        "Calculations" => "calcs.md",
        "API" => "api.md",
        ]
    )
