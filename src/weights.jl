################################################################################
# Goal: compute the average (across instruments) weights on MTRs
#
# The average weights are piecewise-constant functions of u.
# The `compute_average_weights` function will return a data frame that stores
# value of the average weight at the endpoints of each nontrivial piece.
#   1. value of u
#   2. average weight for d = 1
#   3. average weight for d = 0
#
# The average weights for IV-like estimands rely on the s(d, z) function.
# The average weights for target parameters are derived manually.
#
################################################################################

function compute_average_weights(tp::TargetParameter)
    if tp.name == "LATE(u₁, u₂)"
        u₁ = tp.int_limits(1)[1] # only coded the case of 1 instrument
        u₂ = tp.int_limits(1)[2] # only coded the case of 1 instrument
        results = DataFrame(u = [u₁, u₂])
        results[:, "average weight for d = 1"] .= 1 / (u₂ - u₁)
        results[:, "average weight for d = 0"] = -results[:, 2]
    else
        print("Weight Computation is unsupported for this target parameter.")
    end
    return results
end
export compute_average_weights

function compute_average_weights(ivlike::IVLike, dgp::DGP)
    results = DataFrame(u = dgp.pscore)
    order = sortperm(dgp.pscore)
    s = ivlike.s[1] # only coded the case of 1 instrument
    terms = d -> [s(d, dgp.suppZ[i]) for i in 1:length(order)] .* dgp.densZ
    # d = 0
    d0terms = terms(0)
    summands = d0terms[order] # order by increasing pscore
    results[:, "average weight for d = 0"] = cumsum(summands)
    # d = 1
    d1terms = terms(1)
    summands = d1terms[reverse(order)] # order by decreasing pscore
    results[:, "average weight for d = 1"] = reverse(cumsum(summands))
    return results
end
