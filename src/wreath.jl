include("permutations.jl")


# ======== Wreath Products ======== #


"""
    Wreath{N, T}

Wreath product `T ≀ G ≅ Tᴺ ⋊ G`, where `G` is either `Σ{N}` or `Cyc{N}`.
"""
struct Wreath{N, T}
    tuple::NTuple{N, T}
    perm::Perm{N}

    function Wreath{N, T}(tuple::NTuple{N, T}, perm::Perm{N}) where {N, T}
        new{N, T}(tuple, perm)
    end

    function Wreath(tuple::NTuple{N, T}, perm::Perm{N}) where {N, T}
        new{N, T}(tuple, perm)
    end
end

Base.one(::Wreath{N, T}) where {N, T} = Wreath{N, T}(Tuple(one(T) for i in 1:N), one(Σ{N}))
Base.one(::Type{Wreath{N, T}}) where {N, T} = Wreath{N, T}(Tuple(one(T) for i in 1:N), one(Σ{N}))

Base.:(==)(w1::Wreath{N, T}, w2::Wreath{N, T}) where {N, T} = w1.tuple == w2.tuple && w1.perm == w2.perm

function Base.show(io::IO, w::Wreath{N, T}) where {N, T}
    w.tuple == Tuple(one(T) for i in 1:N) ? "" : print(io, w.tuple)
    print(io, w.perm)
end

"""
    ((tuple::NTuple{N, T})^(perm::Perm{N}))::NTuple{N, T}

Permutes `tuple` by `perm`. The resulting tuple satisfies `(tuple^perm)[i] == tuple[i^perm]`.
"""
function Base.:^(tuple::NTuple{N, T}, perm::Perm{N})::NTuple{N, T} where {N, T}
    Tuple(tuple[i^perm] for i in 1:N)
end

"""
    ((t1::NTuple{N, T})*(t2::NTuple{N, T}))::NTuple{N, T}

Pointwise multiplication of tuples.
"""
function Base.:*(t1::NTuple{N, T}, t2::NTuple{N, T}) where {N, T}
    Tuple(t1[i] * t2[i] for i in 1:N)
end

"""
    Base.inv(t::NTuple{N, T})::NTuple{N, T}

Pointwise inverse of tuples.
"""
function Base.inv(t::NTuple{N, T})::NTuple{N, T} where {N, T}
    Tuple(Base.inv(t[i]) for i in 1:N)
end

"""
    ((t::NTuple{N, T})*(n::Int))::NTuple{N, T}

Raises the tuple `t` pointwise to the `n`-th power..
"""
function Base.:^(t::NTuple{N, T}, n::Int)::NTuple{N, T} where {N, T}
    Tuple(t[i]^n for i in 1:N)
end

"""
    ((w1::Wreath{N, T})*(w2::Wreath{N, T}))::Wreath{N, T}

Multiplication in the wreath product.
"""
function Base.:*(w1::Wreath{N, T}, w2::Wreath{N, T})::Wreath{N, T} where {N, T}
    Wreath{N, T}(w1.tuple * (w2.tuple ^ (inv(w1.perm))), w1.perm * w2.perm)
end

"""
    Base.inv(w::Wreath{N, T})::Wreath{N, T}

Taking inverses in the wreath product.
"""
function Base.inv(w::Wreath{N, T})::Wreath{N, T} where {N, T}
    Wreath{N, T}((inv(w.tuple))^(w.perm), inv(w.perm))
end


# ======== Iterated Wreath Products ======== #


"""
    IteratedWreath{N}

Represents an iterated wreath product: `IteratedWreath{N} ≅ IteratedWreath{N} ≀ Σ{N}`.
If `value` is set to `nothing`, it is the identity element.
"""
struct IteratedWreath{N}
    value::Union{Wreath{N, IteratedWreath{N}}, Nothing}

    function IteratedWreath{N}(value::Union{Wreath{N, IteratedWreath{N}}, Nothing}) where N
        new{N}(value)
    end

    function IteratedWreath(tuple::NTuple{N, IteratedWreath{N}}, g::Perm{N}) where N
        new{N}(Wreath{N, IteratedWreath{N}}(tuple, g))
    end
end

Base.one(::IteratedWreath{N}) where N = IteratedWreath{N}(nothing)
Base.one(::Type{IteratedWreath{N}}) where N = IteratedWreath{N}(nothing)

function Base.show(io::IO, w::IteratedWreath{N}) where N
    w.value == nothing ? print(io, "1") : print(io, w.value)
end

function Base.:*(w1::IteratedWreath{N}, w2::IteratedWreath{N}) where N
    w1.value == nothing && return w2
    w2.value == nothing && return w1

    IteratedWreath{N}(w1.value * w2.value)
end

function Base.inv(w::IteratedWreath{N}) where N
    w.value == nothing && return w
    IteratedWreath{N}(inv(w.value))
end
