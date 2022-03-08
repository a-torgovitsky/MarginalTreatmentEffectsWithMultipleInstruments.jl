function menu(savelocation::String = "."; compile::Bool = false)
    println('='^80)
    println("Policy Evaluation with Multiple Instrumental Variables")
    println("Mogstad, Torgovitsky, and Walters (2021)")
    println('='^80)
    println("What do you want to run?")
    println("\t 1. Illustration of mutual consistency (Figure 2)")
    println("\t 2. Numerical illustration for ATT (Figure 4)")
    println("\t 3. Numerical illustration for LATE(+20) (Figure 5)")
    println("\t 4. PRTE misspecification example (Figure 6)")
    println("\t 5. Everything")
    print("Enter choice: ")
    choice = readline()
    choice = parse(Int64, choice)

    if choice == 1
        savedir, _ = setup(savelocation, stub = "illustrate-mc")
        run_illustrate_mc(savedir, compile)
    elseif choice == 2
        savedir, _ = setup(savelocation, stub = "simulation-att")
        run_simulation_att(savedir, compile)
    elseif choice == 3
        savedir, _ = setup(savelocation, stub = "simulation-prte")
        run_simulation_prte(savedir, compile)
    elseif choice == 4
        savedir, _ = setup(savelocation, stub = "prte-misspecification")
        run_prte_misspecification(savedir, compile)
    elseif choice == 5
        savedir, _ = setup(savelocation, stub = "everything")
        results_mc = run_illustrate_mc(savedir, compile)
        results_att = run_simulation_att(savedir, compile)
        results_prte = run_simulation_prte(savedir, compile)
        results_prte_misspecification =
            run_prte_misspecification(savedir, compile)
        return Dict(:results_mc => results_mc,
                    :results_att => results_att,
                    :results_prte => results_prte,
                    :results_prte_misspecification =>
                        results_prte_misspecification)
    end
end
export menu

function run_illustrate_mc(savedir::String, compile::Bool = false)
    results = illustrate_mc()
    fn_results = joinpath(savedir, "illustrate-mc.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigureMTRConsistent.tex"))
        compile_latex(joinpath(savedir, "FigureMTRInconsistent.tex"))
    end
    return results
end
export run_illustrate_mc

function run_simulation_att(savedir::String, compile::Bool = false)
    results = simulation_att()
    fn_results = joinpath(savedir, "simulation-att.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigureATT.tex"))
    end
    return results
end
export run_simulation_att

function run_simulation_prte(savedir::String, compile::Bool = false)
    results = simulation_prte()
    fn_results = joinpath(savedir, "simulation-prte.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigurePRTE.tex"))
    end
    return results
end
export run_simulation_prte

function run_prte_misspecification(savedir::String, compile::Bool = false)
    results = prte_misspecification()
    fn_results = joinpath(savedir, "prte-misspecification.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigurePRTEMisspecification.tex"))
    end
    return results
end
export run_simulation_prte

function compile_latex(fn::String)
    oldwd = pwd()
    try
        cd(dirname(fn))
        cstr = `pdflatex -halt-on-error $(basename(fn)) "|" grep -a3 ^!`
        @suppress begin
            run(cstr)
            run(`latexmk -c`)
        end
        cd(oldwd)
    catch err
        cd(oldwd)
    end
end
