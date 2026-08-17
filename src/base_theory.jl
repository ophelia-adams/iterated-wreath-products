using Metatheory, TermInterface
using Metatheory.Rewriters
using Metatheory.EGraphs

include("wreathexpr.jl")


# ======== Base Theory (Rewrites) ======== #


# Helper function until I figure out how to match on parametrized types.
wreathdim(w::WreathExpr{N}) where N = N

function flattenMul(m::MulExpr{N})::Vector{WreathExpr{N}} where N
    [flattenMul(m.arg1); flattenMul(m.arg2)]
end

function flattenMul(e::WreathExpr{N})::Vector{WreathExpr{N}} where N
    [e]
end

function reduceMul(es::Vector{WreathExpr{N}})::Vector{WreathExpr{N}} where N
    pass1 = reduce_consts(es)

    pass2_0 = nothing
    pass2   = pass1
    while pass2 != pass2_0
        pass2_0 = pass2
        pass2 = elim_inverses(pass2)
    end

    return pass2
end

function elim_inverses(es::Vector{WreathExpr{N}})::Vector{WreathExpr{N}} where N
    length(es) <= 1 && return es
    res = []
    i = 1
    while i < length(es)
        if es[i] == inv(es[i+1])
            i += 2
        else
            i += 1
            push!(res, es[i])
        end
    end
    i == length(es) && push!(res, es[i])
    res
end

function reduce_consts(es::Vector{WreathExpr{N}})::Vector{WreathExpr{N}} where N
    length(es) <= 1 && return es
    res = []
    buildup = nothing
    i = 1
    while i <= length(es)
        if es[i] isa(ConExpr)
            if buildup == nothing
                buildup = es[i]
            else
                buildup *= es[i]
            end
        else
            if buildup != nothing
                push!(res, buildup)
                buildup = nothing
            end
            push!(res, es[i])
        end
        i += 1
    end
    res
end

function reconstructMul(es::Vector{WreathExpr{N}})::WreathExpr{N} where N
    isempty(es) && return ConExpr(one(IteratedWreath{N}))
    length(es) == 1 && return es[1]
    *(es...)
end


# phase1: basic groups theory
theory_phase1 = @theory x begin
    x * inv(x) => ConExpr(one(IteratedWreath{wreathdim(x)}))
    inv(x) * x => ConExpr(one(IteratedWreath{wreathdim(x)}))
    inv(inv(x)) => x
end

r_phase1 = Fixpoint(Postwalk(Chain(theory_phase1)))

# phase2: flatten multiplication to clean it up more
r_phase2 = @rule ~x => reconstructMul(reduceMul(flattenMul(~x)))


# ======== Demo ======== #

#=
You need Julia 1.11 for this, as 1.12 seems to break MatchCore,
on which Metatheory and friends depend

In the REPL, include this file, and define the following variable and test expression:

julia> s = VarExpr{3}(:s)
VarExpr{3}(:s)

julia> ex = s * inv(s) * inv(s) * inv(s) * (s * inv(s)) * s * s
MulExpr{3}(MulExpr{3}(MulExpr{3}(MulExpr{3}(MulExpr{3}(MulExpr{3}(VarExpr{3}(:s), InvExpr{3}(VarExpr{3}(:s))), InvExpr{3}(VarExpr{3}(:s))), InvExpr{3}(VarExpr{3}(:s))), MulExpr{3}(VarExpr{3}(:s), InvExpr{3}(VarExpr{3}(:s)))), VarExpr{3}(:s)), VarExpr{3}(:s))

The rule `r_phase1` traverses the expression and its subexpressions,
    looking for multipliation by inverses and inverses of inverses.
However, as many instances of `s` and `inv(s)` are linearly adjacent but structurally far,
   this rule does not get us far.

julia> r_phase1(ex)
MulExpr{3}(MulExpr{3}(MulExpr{3}(InvExpr{3}(VarExpr{3}(:s)), InvExpr{3}(VarExpr{3}(:s))), VarExpr{3}(:s)), VarExpr{3}(:s))

The rule `r_phrase2`, on the other hand, flattens the product completely (see the functions above):

julia> r_phase2(ans)
ConExpr{3}(1)

It is evidently more effective here, but the trade-off is that the walkers provided by Metatheory.Rewriters
    are not efficient here (as they would have to flatten, reduce, unflatten, flatten, reduce, etc.).

=#
