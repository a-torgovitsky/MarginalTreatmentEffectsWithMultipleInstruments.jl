# Mogstad, Walters, Torgovitsky (2021)

Repo: https://github.com/a-torgovitsky/MarginalTreatmentEffectsWithMultipleInstruments.jl

This package contains the code for the figures and simulations in "Policy Evaluation with Multiple Instrumental Variables" by Mogstad, Walters, and Torgovitsky (2021).
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
