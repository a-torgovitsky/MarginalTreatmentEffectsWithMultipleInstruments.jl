function make_slist(suppZ)
    # one way to generalize is allow for more general slists
    # for these simulations we are only using fully saturated
    return [((d,z) -> Int((d == d̄) * (z == z̄)))
            for d̄ in 0:1, z̄ in eachrow(suppZ)][:]
end

function compute_βₛ(dgp::DGP; slist = "saturated")
    Γₛ = compute_Γₛ([(dgp.mtrs[1].basis, dgp.mtrs[2].basis)], dgp,
                    slist = slist)
    # Only one model here, so Γₛ is a 1 x 2 array of matrices
    βₛ = fill(NaN, size(Γₛ[1,1])[1])
    for s in 1:length(βₛ)
        βₛ[s] = sum(Γₛ[1,1][s,:,:] .* dgp.mtrs[1].θ +
                    Γₛ[1,2][s,:,:] .* dgp.mtrs[2].θ)
    end
    return βₛ
end
export compute_βₛ

# return a 2-dimensional array Γₛ[ℓ, d], where each component is itself
# a 3-dimensional array with elements [s, j, k]
function compute_Γₛ(
    bases::Array{Tuple{MTRBasis, MTRBasis}, 1},
    dgp::DGP;
    slist = "saturated"
)
    [compute_Γₛ(basis[d + 1], d, dgp, slist = slist) for basis in bases,
                                                         d in 0:1]
end

function compute_Γₛ(basis::MTRBasis, d::Integer, dgp::DGP; slist = "saturated")
    @assert d in [0,1]
    if (slist == "saturated")
        slist = make_slist(dgp.suppZ)
    end
    Γₛ = zeros(length(slist), length(basis.a), length(basis.b))
    for (i,z) in enumerate(eachrow(dgp.suppZ))
        intlb = (1 - d) * dgp.pscore[i]
        intub = d * dgp.pscore[i] + (1 - d) * 1
        Γₛ += [(aj(z) * ibk(intlb, intub) * s(d, z))
                for s in slist, aj in basis.a, ibk in basis.ib] .* dgp.densZ[i]
    end
    return Γₛ
end
