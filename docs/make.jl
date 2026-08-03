using Documenter
using MyPackage

makedocs(
    sitename = "Aircraft Optimizer",
    modules = [MyPackage],

    format = Documenter.HTML(
        inventory_version = string(pkgversion(MyPackage)),
        assets = ["assets/style.css"],
    ),

    pages = [
        "Home" => "index.md",
        "Getting Started" => "onboarding/getting_started.md",
        "Understanding the Repository" => "onboarding/repo_structure.md",
        "Technical" => [
            "Sign Conventions" => "technical/sign_conventions.md"
        ],
        "Code Reference" => "code_reference/code_reference.md",
    ],
)

deploydocs(
    repo="github.com/ArthurJWH/Aircraft_Optimizer.git",
    devbranch="main",
)