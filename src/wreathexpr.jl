using Metatheory, TermInterface
using Metatheory.EGraphs

include("wreath.jl")


# ======== Expression Types (Syntax Tree Components) ======== #


abstract type WreathExpr{N, T} end

struct ConExpr{N, T} <: WreathExpr{N, T}
    value::Wreath{N, T}
end

struct VarExpr{N, T} <: WreathExpr{N, T}
    name::Symbol
end

struct TupExpr{N, T} <: WreathExpr{N, T}
    tuple::NTuple{N, WreathExpr{N, T}}
end

struct MulExpr{N, T} <: WreathExpr{N, T}
    arg1::WreathExpr{N, T}
    arg2::WreathExpr{N, T}
end

struct InvExpr{N, T} <: WreathExpr{N, T}
    arg::WreathExpr{N, T}
end


# ======== TermInterface Methods ======== #


TermInterface.istree(::WreathExpr) = true
TermInterface.istree(::ConExpr) = false
TermInterface.istree(::VarExpr) = false


function Base.:*(t1::NTuple{N, T1}, t2::NTuple{N, T2}) where {N, T1 <: WreathExpr, T2 <: WreathExpr}
    Tuple(t1[i] * t2[i] for i in 1:N)
end

Base.:*(w1::WreathExpr, w2::WreathExpr) = MulExpr(w1, w2)
Base.:*(c1::ConExpr, c2::ConExpr) = ConExpr(c1.value * c2.value)

function Base.:*(c::ConExpr{N, T}, t::TupExpr{N, T}) where {N, T}
    permuted_tuple = t.tuple^((c.value.perm)^(-1))
    consted_tuple = map(ConExpr{N, T}, c.value.tuple)
    res_tuple = consted_tuple * permuted_tuple

    res_perm = c.value.perm
    res_const = Wreath{N, T}(Tuple(res_perm^0 for i in 1:N), res_perm)

    MulExpr{N, T}(TupExpr{N, T}(res_tuple), ConExpr{N, T}(res_const))
end

Base.inv(w::WreathExpr) = InvExpr(w)
Base.inv(c::ConExpr) = ConExpr(inv(c.value))

TermInterface.operation(e::TupExpr) = tuple
TermInterface.operation(e::MulExpr) = *
TermInterface.operation(e::InvExpr) = inv


TermInterface.arguments(e::TupExpr) = collect(e.tuple)
TermInterface.arguments(e::MulExpr) = [e.arg1, e.arg2]
TermInterface.arguments(e::InvExpr) = [e.arg]


TermInterface.exprhead(::WreathExpr) = :call

TermInterface.metadata(::WreathExpr) = nothing


#=
function TermInterface.similarterm(x::ConExpr, head, args; metadata = nothing, exprhead = :call)
    ConExpr(args[1])
end

function TermInterface.similarterm(x::VarExpr, head, args; metadata = nothing, exprhead = :call)
    VarExpr(args[1])
end
=#

function TermInterface.similarterm(x::TupExpr, head, args; metadata = nothing, exprhead = :call)
    TupExpr(tupleop(args...))
end

function TermInterface.similarterm(x::MulExpr, head, args; metadata = nothing, exprhead = :call)
    MulExpr(args[1], args[2])
end

function TermInterface.similarterm(x::InvExpr, head, args; metadata = nothing, exprhead = :call)
    InvExpr(args[1])
end


#=
function EGraphs.egraph_reconstruct_expression(::Type{ConExpr{N, T}}, op, args; metadata = nothing, exprhead = nothing) where {N, T}
    ConExpr{N, T}(args[1])
end

function EGraphs.egraph_reconstruct_expression(::Type{VarExpr{N, T}}, op, args; metadata = nothing, exprhead = nothing) where {N, T}
    VarExpr{N, T}(args[1])
end
=#

function EGraphs.egraph_reconstruct_expression(::Type{TupExpr{N, T}}, op, args; metadata = nothing, exprhead = nothing) where {N, T}
    TupExpr{N, T}(tupleop(args...))
end

function EGraphs.egraph_reconstruct_expression(::Type{MulExpr{N, T}}, op, args; metadata = nothing, exprhead = nothing) where {N, T}
    MulExpr{N, T}(args[1], args[2])
end

function EGraphs.egraph_reconstruct_expression(::Type{InvExpr{N, T}}, op, args; metadata = nothing, exprhead = nothing) where {N, T}
    InvExpr{N, T}(args[1])
end


# ======== Base Theory (Rewrites) ======== #


# Helper functions until I figure out how to match on parametrized types.
wreathdim(w::WreathExpr{N, T}) where {N, T} = N
wreathtype(w::WreathExpr{N, T}) where {N, T} = T

theory_groups = @theory x begin
    x * inv(x) => VarExpr{wreathdim(x), wreathtype(x)}(:id)
    inv(x) * x => VarExpr{wreathdim(x), wreathtype(x)}(:id)
    inv(inv(x)) => x
    
end
