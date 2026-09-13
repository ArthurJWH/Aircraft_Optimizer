using Documenter
using MyPackage
using MyPackage.Utils
using MyPackage.PlaneInfo
using MyPackage.Geometry
using MyPackage.MyPlots
using MyPackage.VLM
using MyPackage.IO

makedocs(
    sitename = "Aircraft Optimizer",
    authors = "Arthur JWH",
    modules = [
        MyPackage,
        MyPackage.Utils,
        MyPackage.PlaneInfo,
        MyPackage.Geometry,
        MyPackage.MyPlots,
        MyPackage.VLM,
        MyPackage.IO,
    ],

    format = Documenter.HTML(
        inventory_version = string(pkgversion(MyPackage)),
        assets = [
            "assets/style.css",
            "assets/favicon.ico",
        ],
        description = "Aircraft Optimizer Documentation — A VLM based workflow"
    ),

    pages = [
        "Home" => "index.md",
        "Getting Started" => "onboarding/getting_started.md",
        "Understanding the Repository" => "onboarding/repo_structure.md",
        "User Guide" => [
            "Your First Simulation" => "how_to/first_simulation.md",
            "Post-Processing" => "how_to/post_processing.md",
            "Stability Analysis" => "how_to/stability_analysis.md",
            "Ground Effect" => "how_to/ground_effect.md",
            "CAD Export" => "how_to/cad_export.md",
        ],
        "Technical" => [
            "Sign Conventions" => "technical/sign_conventions.md",
            "VLM Theory" => "technical/vlm_theory.md",
            "Solver Architecture" => "technical/solver_architecture.md",
        ],
        "Code Reference" => "code_reference/code_reference.md",
    ],
)

deploydocs(
    repo="github.com/ArthurJWH/Aircraft_Optimizer.git",
    devbranch="main",
)