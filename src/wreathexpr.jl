using Metatheory, TermInterface
using Metatheory.EGraphs

include("wreath.jl")


# ======== Expression Types (Syntax Tree Components) ======== #


abstract type WreathExpr{N} end

struct ConExpr{N} <: WreathExpr{N}
    value::IteratedWreath{N}

    ConExpr(w::IteratedWreath{N}) where N = new{N}(w)
end

struct VarExpr{N} <: WreathExpr{N}
    name::Symbol
end

struct TupExpr{N} <: WreathExpr{N}
    tuple::NTuple{N, WreathExpr{N}}

    TupExpr(tuple::NTuple{N}) where N = new{N}(tuple)
end

struct MulExpr{N} <: WreathExpr{N}
    arg1::WreathExpr{N}
    arg2::WreathExpr{N}

    MulExpr(arg1::WreathExpr{N}, arg2::WreathExpr{N}) where N = new{N}(arg1, arg2)
end

struct InvExpr{N} <: WreathExpr{N}
    arg::WreathExpr{N}

    InvExpr(w::WreathExpr{N}) where N = new{N}(w)
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

function Base.:*(c::ConExpr{N}, t::TupExpr{N}) where N
    c.value.value == nothing && return t

    permuted_tuple = t.tuple^((c.value.value.perm)^(-1))
    consted_tuple = map(ConExpr{N}, c.value.value.tuple)
    res_tuple = consted_tuple * permuted_tuple

    res_perm = c.value.value.perm
    res_const = IteratedWreath{N}(Tuple(res_perm^0 for i in 1:N), res_perm)

    MulExpr{N}(TupExpr{N}(res_tuple), ConExpr{N}(res_const))
end

function Base.:*(c::ConExpr{N}, w::WreathExpr{N}) where N
    c.value.value == nothing ? w : MulExpr(c, w)
end

function Base.:*(w::WreathExpr{N}, c::ConExpr{N}) where N
    c.value.value == nothing ? w : MulExpr(w, c)
end

Base.inv(w::WreathExpr) = InvExpr(w)
Base.inv(c::ConExpr) = ConExpr(inv(c.value))
Base.inv(w::InvExpr) = w.arg

TermInterface.operation(e::TupExpr) = :tuple
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
    TupExpr(tuple(args...))
end

function TermInterface.similarterm(x::MulExpr, head, args; metadata = nothing, exprhead = :call)
    args[1] * args[2]
end

function TermInterface.similarterm(x::InvExpr, head, args; metadata = nothing, exprhead = :call)
    inv(args[1])
end


#=
function EGraphs.egraph_reconstruct_expression(::Type{ConExpr{N}}, op, args; metadata = nothing, exprhead = nothing) where N
    ConExpr{N}(args[1])
end

function EGraphs.egraph_reconstruct_expression(::Type{VarExpr{N}}, op, args; metadata = nothing, exprhead = nothing) where N
    VarExpr{N}(args[1])
end
=#

function EGraphs.egraph_reconstruct_expression(::Type{TupExpr{N}}, op, args; metadata = nothing, exprhead = nothing) where N
    TupExpr{N}(tupleop(args...))
end

function EGraphs.egraph_reconstruct_expression(::Type{MulExpr{N}}, op, args; metadata = nothing, exprhead = nothing) where N
    args[1] * args[2]
end

function EGraphs.egraph_reconstruct_expression(::Type{InvExpr{N}}, op, args; metadata = nothing, exprhead = nothing) where N
    inv(args[1])
end


