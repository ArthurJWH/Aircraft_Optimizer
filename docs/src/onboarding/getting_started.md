# Getting Started

## Set up

1. Download [Visual Studio Code](https://code.visualstudio.com/download?_exp_download=fb315fc982) (recommended source code editor)

2. Download [Git](https://git-scm.com/install/) for version control.

3. Download [Julia](https://julialang.org/downloads/manual-downloads/) ("Long-term Support Release" recommended): a Just-in-time (JIT) compilation programming language with syntax similar to Python and Matlab.

    !!! warning "Recall"
        During installation, in "Select Additional Tasks", check "Add Julia to PATH".

4. Open Visual Studio Code, go to Terminal > New Terminal, and run the following lines:

```
git config --global user.name "Somebody"
git config --global user.email somebody@mail.com
```

5. Open EXTENSIONS tab and install:
    - Julia (identifier: julialang.language-julia)
    - GitHub Repository Manager (identifier: henriquebruno.github-repository-manager)

6. Open GitHub Repository Manager tab:
    - Log in into your GitHub account.
    - Click on Aircraft Optimizer repository and select a folder to clone the source code (recommended: "C:\").

7. Go to File > Open Folder..., and select the cloned folder.

8. Open Setup.jl and run it once by pressing Alt+Enter.

!!! note "Note"
    At the start of each session, it might be required to run Startup.jl