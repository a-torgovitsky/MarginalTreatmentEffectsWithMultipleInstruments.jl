# Illustrate the concept of mutual consistency
function illustrate_mc()
    dgp = illustration_dgp()
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1)
    basis1 = (interacted_bernstein_basis(1, notℓ = 2),
              interacted_bernstein_basis(1, notℓ = 2))
    basis2 = (interacted_bernstein_basis(3, notℓ = 1),
              interacted_bernstein_basis(3, notℓ = 1))
    bases = [basis1, basis2]
    tp = ey1(dgp, ℓ = 2)

    # for plotting
    ev = DataFrame(u = 0:.01:1, d = 1)
    results = DataFrame(u = ev.u)

    # Impose mutual consistency
    assumptions[:mutually_consistent] = true
    r = compute_bounds(tp, bases, assumptions, dgp)
    ev.z = fill([0,0], nrow(ev))
    results[:, "MTR (l=1 d=1 z2=0)"] = evaluate_mtr(r[:mtr_ub][1][2], ev)
    ev.z = fill([0,1], nrow(ev))
    results[:, "MTR (l=1 d=1 z2=1)"] = evaluate_mtr(r[:mtr_ub][1][2], ev)
    ev.z = fill([0,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=0) con"] = evaluate_mtr(r[:mtr_ub][2][2], ev)
    ev.z = fill([1,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=1) con"] = evaluate_mtr(r[:mtr_ub][2][2], ev)

    # Do not impose mutual consistency. Keep target parameter defined in terms
    # of the second model. Then the max MTR will be a second model MTR that is
    # inconsistent with the first model MTR.
    assumptions[:mutually_consistent] = false
    r = compute_bounds(tp, bases, assumptions, dgp)
    ev.z = fill([0,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=0) incon"] = evaluate_mtr(r[:mtr_ub][2][2], ev)
    ev.z = fill([1,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=1) incon"] = evaluate_mtr(r[:mtr_ub][2][2], ev)
    return results
end
export illustrate_mc

# Show the impact of mutual consistency for an instrument-invariant parameter
# like the ATT
function simulation_att()
    dgp = simulation_dgp()

    results = DataFrame(degree = 1:1:9)
    results[:, "LB (l = 1; MC off)"] .= NaN
    results[:, "UB (l = 1; MC off)"] .= NaN
    results[:, "LB (l = 2; MC off)"] .= NaN
    results[:, "UB (l = 2; MC off)"] .= NaN
    results[:, "LB (MC on)"] .= NaN
    results[:, "UB (MC on)"] .= NaN
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1)

    for k in 1:nrow(results)
        bases = [(interacted_bernstein_basis(results[k, :degree], notℓ = 2),
                  interacted_bernstein_basis(results[k, :degree], notℓ = 2)),
                 (interacted_bernstein_basis(results[k, :degree], notℓ = 1),
                  interacted_bernstein_basis(results[k, :degree], notℓ = 1))]

        assumptions[:mutually_consistent] = false
        tp = att(dgp, ℓ = 1)
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 2:3] = [r[:lb], r[:ub]]

        tp = att(dgp, ℓ = 2)
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 4:5] = [r[:lb], r[:ub]]

        assumptions[:mutually_consistent] = true
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 6:7] = [r[:lb], r[:ub]]
    end

    knots = vcat(0, dgp.pscore, 1)
    bases = [(interacted_constantspline_basis(knots, notℓ = 2),
              interacted_constantspline_basis(knots, notℓ = 2)),
             (interacted_constantspline_basis(knots, notℓ = 1),
              interacted_constantspline_basis(knots, notℓ = 1))]

    assumptions[:mutually_consistent] = false
    tp = att(dgp, ℓ = 1)
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[:, "NP LB (l = 1; MC off)"] .= r[:lb]
    results[:, "NP UB (l = 1; MC off)"] .= r[:ub]

    tp = att(dgp, ℓ = 2)
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[:, "NP LB (l = 2; MC off)"] .= r[:lb]
    results[:, "NP UB (l = 2; MC off)"] .= r[:ub]

    assumptions[:mutually_consistent] = true
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[:, "NP LB (MC on)"] .= r[:lb]
    results[:, "NP UB (MC on)"] .= r[:ub]

    return results
end
export simulation_att

# Show that mutual consistency can also help tighten inference on instrument-
# dependent parameters.
function simulation_prte()
    dgp = simulation_dgp()
    δ = .2
    tp = prte_plusδpercent(dgp, δ, ℓ = 1)
    knots = vcat(0, dgp.pscore, (1 + δ) .* dgp.pscore, 1)

    results = DataFrame(degree = 1:1:9)
    results[:, "LB (MC off)"] .= NaN
    results[:, "UB (MC off)"] .= NaN
    results[:, "LB (l = 2; NP)"] .= NaN
    results[:, "UB (l = 2; NP)"] .= NaN
    results[:, "LB (l = 2; NP; Decr)"] .= NaN
    results[:, "UB (l = 2; NP; Decr)"] .= NaN
    results[:, "LB (l = 2; NP; Linear)"] .= NaN
    results[:, "UB (l = 2; NP; Linear)"] .= NaN

    for k in 1:nrow(results)
        bases = [(interacted_bernstein_basis(results[k, :degree], notℓ = 2),
                  interacted_bernstein_basis(results[k, :degree], notℓ = 2)),
                 (interacted_constantspline_basis(knots, notℓ = 1),
                  interacted_constantspline_basis(knots, notℓ = 1))]

        assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1,
                                        :mutually_consistent => false)
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 2:3] = [r[:lb], r[:ub]]

        assumptions[:mutually_consistent] = true
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 4:5] = [r[:lb], r[:ub]]

        assumptions[:decreasing_level] = [(2, 0), (2,1)]
        assumptions[:decreasing_difference] = [2]
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 6:7] = [r[:lb], r[:ub]]

        bases[2] = (interacted_bernstein_basis(1, notℓ = 1),
                    interacted_bernstein_basis(1, notℓ = 1))
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 8:9] = [r[:lb], r[:ub]]
    end

    return results
end
export simulation_prte

# Illustrate misspecification caused for a PRTE by using the wrong choice model
function prte_misspecification(; ℓ_gen = 1)
    dgp = prte_dgp(ℓ_gen = ℓ_gen) # ℓ_gen shouldn't matter here (and doesn't)
    knots = vcat(0, dgp.pscore, 1)

    # Set up PRTE target parameter
    dgp_new = DGP(suppZ = dgp.suppZ,
                  pscore = dgp.pscore,
                  densZ = fill(.25, 4),
                  mtrs = dgp.mtrs)
    tp = prte_newz(dgp, dgp_new, ℓ = 1)

    nrows = 5
    results = DataFrame(name = fill("", nrows),
                        lb = fill(+Inf, nrows), ub = fill(+Inf, nrows))

    # True value
    truth = eval_tp(tp, [dgp.mtrs], dgp)
    results[1, :name] = "True value"
    results[1, :lb], results[1, :ub] = truth, truth

    # IAM model using both instruments
    # --> should get bounds that are non-positive
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1)
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[2, :name] = "Both instruments with IAM"
    results[2, :lb] = r[:lb]
    results[2, :ub] = r[:ub]

    # Use the first instrument only
    bases = [(interacted_constantspline_basis(knots, notℓ = 2),
              interacted_constantspline_basis(knots, notℓ = 2))]
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[3, :name] = "Instrument 1 with IAM"
    results[3, :lb] = r[:lb]
    results[3, :ub] = r[:ub]

    # Use the second instrument only
    bases = [(interacted_constantspline_basis(knots, notℓ = 1),
              interacted_constantspline_basis(knots, notℓ = 1))]
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[4, :name] = "Instrument 2 with IAM"
    results[4, :lb] = r[:lb]
    results[4, :ub] = r[:ub]

    # Use both instruments together and impose MC
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1,
                                    :mutually_consistent => true)
    bases = [(interacted_constantspline_basis(knots, notℓ = 2),
              interacted_constantspline_basis(knots, notℓ = 2)),
             (interacted_constantspline_basis(knots, notℓ = 1),
              interacted_constantspline_basis(knots, notℓ = 1))]
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[5, :name] = "Both instruments with PM"
    results[5, :lb] = r[:lb]
    results[5, :ub] = r[:ub]
    return results
end
export prte_misspecification
