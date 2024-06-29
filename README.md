# MMPSingularities.jl

This package aims to implement cutting-edge-speed algorithms that calculate various invariants of singularities, especially those singularities that show up in the Minimal Model Program. 

**This Package is a work in progress, under active development.** Contributions and PRs are welcome. If you want to add an algorithm that hasn't been implemented, please reach out via email to [Jack J. Garzella](https://mathweb.ucsd.edu/~jjgarzel/) to make sure we don't duplicate work.

This package would not be possible without CUDA.jl, and the broader Julia GPU community in general. Moreover, this package is inspired by the Macaulay2 computational algebra system and especially the community behind it. 

## Current functionality

### Quasi-F-Split Height

We currently implement cpu-only and gpu-accelerated versions of

* Fedder's criterion for Calabi-Yau hypersurfaces (see Theorem C of [arXiv:2204.10076](https://arxiv.org/abs/2204.10076))

```
TODO: sample code here
```

We have cpu-only versions of

* Fedder's criterion for prinicpal ideals (see Theorem A of [arXiv:2204.10076](https://arxiv.org/abs/2204.10076)
* The classical Fedder's criterion for principal ideals (see [Fedder's original paper](https://www.ams.org/journals/tran/1983-278-02/S0002-9947-1983-0701505-0/S0002-9947-1983-0701505-0.pdf))

```
TODO: sample code here
```

### Thresholds

We implement

* the LCT/FPT of determinantal varieties (see e.g. [arXiv:1210.6729](https://arxiv.org/abs/1210.6729))

## Getting Started

Since some singularities researchers might be new to Julia, we don't assume that you already have Julia installed on your system.

1. [Install `juliaup`](https://github.com/JuliaLang/juliaup) using the command line
2. Install your desired Julia version using `juliaup`
3. (Optional, but recommended) Install [VS Code](https://code.visualstudio.com/) and the [Julia plugin](https://www.julia-vscode.org/)
4. If you have an Nvidia GPU, install the Nvidia driver and CUDA.jl (see the [CUDA.jl docs](https://cuda.juliagpu.org/stable/installation/overview/))
5. Install MMPSingularities.jl:
  * Open a Julia REPL
  * type `]` to enter package mode
  * `add [insert-respository-clone-url-here]`
6. Press backspace to exit package mode, then `using MMPSingularities` (TODO: this won't work right now, fix that)

### Sample code

```
TODO: sample code here
```
