################################################################################
# An IVLike contains a vector of IV-like specifications, one for each
# instrument.
# Each IV-like specification is an anonymous function of two variables:
#   1. treatment status d, and
#   2. corresponding instrument value.
################################################################################
struct IVLike
    name::String
    s::Vector

    function IVLike(name, s)
        new(name, s)
    end
end


function make_slist(suppZ)
    # one way to generalize is allow for more general slists
    # for these simulations we are only using fully saturated
    name = "Saturated"
    s = [((d,z) -> (d == d̄) * (z == z̄))
            for d̄ in 0:1, z̄ in eachrow(suppZ)][:]
    IVLike(name, s)
end

function ivslope(dgp::DGP)
    @assert size(dgp.suppZ, 2) == 1 # haven't coded other cases
    name = "IV Slope"
    expZ = dot(dgp.suppZ, dgp.densZ)
    expD = dot(dgp.pscore, dgp.densZ)
    expDZ = dot(dgp.pscore, dgp.densZ .* dgp.suppZ)
    covDZ = expDZ - expD * expZ
    s = [((d,z) -> ((z[1] - expZ) / covDZ))][:]
    IVLike(name, s)
end
export ivslope

function olsslope(dgp::DGP)
    @assert size(dgp.suppZ, 2) == 1 # haven't coded other cases
    name = "OLS Slope"
    prd1 = dot(dgp.pscore, dgp.densZ)
    s = [((d,z) -> ((d - prd1) / (prd1 * (1 - prd1))))][:]
    IVLike(name, s)
end
export olsslope

function ivslope_indicator(dgp::DGP; support = [1])
    @assert size(dgp.suppZ, 2) == 1 # haven't coded other cases
    set = replace(string(support), "[" => "{", "]" => "}")
    name = "IV Slope for 𝟙(Z == z) for z ∈" * set
    expZind = i -> dgp.densZ[i]
    expDZind = i -> dgp.pscore[i] * dgp.densZ[i]
    expD = dot(dgp.pscore, dgp.densZ)
    covDZind = i -> expDZind(i) - expD * expZind(i)
    s = [((d,z) -> (((z[1] == dgp.suppZ[i]) - expZind(i)) / covDZind(i)))
         for i in support][:]
    IVLike(name, s)
end
export ivslope_indicator

function tslsslope_indicator(dgp::DGP)
    @assert size(dgp.suppZ, 2) == 1 # haven't coded other cases
    set = replace(string(support), "[" => "{", "]" => "}")
    name = "TSLS Slope for 𝟙(Z == z) for z ∈" * set
    Ztilde = dgp.densZ
    expZZ = diagm(Ztilde)
    expDZind = dgp.densZ .* dgp.pscore
    expZX = hcat(Ztilde, expDZind)
    Π = expZX' * inv(expZZ)
    mult = inv(Π * expZX) * Π
    s = [((d,z) -> (mult * [z[1] == zi for zi in dgp.suppZ])[2])][:]
    IVLike(name, s)
end
export tslsslope_indicator

function wald(dgp::DGP; z₀, z₁)
    @assert size(dgp.suppZ, 2) == 1 # haven't coded other cases
    name = "Wald (" * string(z₀) * ", " * string(z₁) * ")"
    dens1 = find_density([z₁], dgp)
    dens0 = find_density([z₀], dgp)
    pscore1 = find_pscore([z₁], dgp)
    pscore0 = find_pscore([z₀], dgp)
    denom = pscore1 - pscore0
    s = [((d,z) -> (((z[1] == z₁)/dens1 - (z[1] == z₀)/dens0) / denom))][:]
    IVLike(name, s)
end
export wald

function compute_βₛ(dgp::DGP; slist = "saturated", param = missing)
    Γₛ = compute_Γₛ([(dgp.mtrs[1].basis, dgp.mtrs[2].basis)], dgp,
                    slist = slist, param = param)
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
    slist = "saturated",
    param = missing
)
    [compute_Γₛ(basis[d + 1], d, dgp, slist = slist, param = param)
     for basis in bases, d in 0:1]
end

function compute_Γₛ(
    basis::MTRBasis,
    d::Integer,
    dgp::DGP;
    slist = "saturated",
    param = missing
)
    @assert d in [0,1]
    if (slist == "saturated")
        slist = make_slist(dgp.suppZ).s
    elseif (slist == "ivslope")
        slist = ivslope(dgp).s
    elseif (slist == "olsslope")
        slist = olsslope(dgp).s
    elseif (slist == "ivslopeind")
        slist = ivslope_indicator(dgp, support = param).s
    elseif (slist == "tslsslopeind")
        slist = tslsslope_indicator(dgp).s
    elseif (slist == "wald")
        slist = []
        for p in param
            push!(slist, wald(dgp; z₀ = p[1], z₁ = p[2]).s[1])
        end
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
