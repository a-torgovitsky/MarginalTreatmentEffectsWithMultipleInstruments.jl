***This package is superceded by [MarginalTreatmentEffects.jl](https://github.com/a-torgovitsky/MarginalTreatmentEffects.jl).***
MarginalTreatmentEffects.jl replicates the figures in the following papers:

1. "Using Instrumental Variables for Inference About Policy Relevant Treatment Parameters"
    Mogstad, Santos, and Torgovitsky (2018, _Econometrica_)
2. "Identification and Extrapolation of Causal Effects with Instrumental Variables"
    Mogstad and Torgovitsky (2018, _Annual Review of Economics_)
3. "Policy Evaluation with Multiple Instrumental Variables"
    Mostad, Torgovitsky, and Walters (Forthcoming, _Journal of Econometrics_)

# Mogstad, Torgovitsky, and Walters (Forthcoming, _Journal of Econometrics_)

Repo: https://github.com/a-torgovitsky/MarginalTreatmentEffectsWithMultipleInstruments.jl

This package contains the code for the figures and simulations in "Policy Evaluation with Multiple Instrumental Variables" by Mogstad, Torgovitsky, and Walters (Forthcoming).
The code is written in Julia and tested with Julia version 1.6.0.
To load/install dependencies:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate() # download the correct packages
```
Then run with
```julia
using MarginalTreatmentEffectsWithMultipleInstruments
menu("/path/to/where/you/want/to/save")
```
If you want to try to automatically call `latexpdf` to compile the TikZ figures, do this:
```julia
menu("/path/to/where/you/want/to/save", compile = true)
```
