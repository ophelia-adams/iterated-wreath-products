# Various conversion, cross-multiplication, etc. functions for ease of use.
# ================================

include("wreathexpr.jl")


Base.convert(::Type{Cyc{N}}, k::Int) where N = Cyc{N}(k)

function Base.convert(::Type{Σ{N}}, g::Cyc{N}) where N
    ncycle = Σ{N}(Tuple(mod1(i+1, N) for i in 1:N), g.name)
    ncycle^g.pow
end

Base.:*(a::Cyc{N}, b::Σ{N}) where N = convert(Σ{N}, a) * b
Base.:*(a::Σ{N}, b::Cyc{N}) where N = a * convert(Σ{N}, b)

Base.:(==)(a::Cyc{N}, b::Σ{N}) where N = convert(Σ{N}, a) == b
Base.:(==)(a::Σ{N}, b::Cyc{N}) where N = a == convert(Σ{N}, b)


function Base.convert(::Type{Wreath{N, T}}, g::Cyc{N}) where {N, T}
    trivial_tuple = Tuple(one(T) for i in 1:N)
    Wreath{N, T}(trivial_tuple, g)
end

function Base.convert(::Type{Wreath{N, T}}, g::Σ{N}) where {N, T}
    trivial_tuple = Tuple(one(T) for i in 1:N)
    Wreath{N, T}(trivial_tuple, g)
end


# For some reason, `g::Perm{N}` did not work at all so I had to split it.
function Base.convert(::Type{IteratedWreath{N}}, g::Cyc{N}) where N
    tuple = NTuple{N}(IteratedWreath{N}(nothing) for i in 1:N)
    IteratedWreath(tuple, g)
end

function Base.convert(::Type{IteratedWreath{N}}, g::Σ{N}) where N
    tuple = NTuple{N}(IteratedWreath{N}(nothing) for i in 1:N)
    IteratedWreath(tuple, g)
end

function Base.convert(::Type{IteratedWreath{N}}, w::Wreath{N, T}) where {N, T}
    convertedtuple = map(x -> x == one(T) ? IteratedWreath{N}(nothing)
                         : Base.convert(IteratedWreath{N}, x), w.tuple)
    IteratedWreath{N}(Wreath{N, IteratedWreath{N}}(convertedtuple, w.perm))
end


Base.convert(::Type{ConExpr{N}}, w::IteratedWreath{N}) where N = ConExpr(w)

function Base.convert(::Type{ConExpr{N}}, w::Wreath{N, T}) where {N, T}
    ConExpr(convert(IteratedWreath{N}, w))
end

