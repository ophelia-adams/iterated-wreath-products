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
end

Base.one(W::Wreath{N,T}) where {N, T} = Wreath{N,T}(Tuple(one(T) for _ in 1:N), one(W.perm))

function Base.show(io::IO, w::Wreath{N, T}) where {N, T}
    print(io, w.tuple)
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
    Tuple(Base.inv(t1[i]) for i in 1:N)
end

"""
    ((t::NTuple{N, T})*(n::Int))::NTuple{N, T}

Raises the tuple `t` pointwise to the `n`-th power..
"""
function Base.:^(t::NTuple{N, T}, n::Int)::NTuple{N, T} where {N, T}
    Tuple(t1[i]^n for i in 1:N)
end

"""
    ((w1::Wreath{N, T})*(w2::Wreath{N, T}))::Wreath{N, T}

Multiplication in the wreath product.
"""
function Base.:*(w1::Wreath{N, T}, w2::Wreath{N, T})::Wreath{N, T} where {N, T}
    Wreath{N, T}(w1.tuple * (w2.tuple ^ (w1.perm^(-1))), w1.perm * w2.perm)
end

"""
    Base.inv(w::Wreath{N, T})::Wreath{N, T}

Taking inverses in the wreath product.
"""
function Base.inv(w::Wreath{N, T})::Wreath{N, T} where {N, T}
    Wreath{N, T}((w.tuple^(-1))^(w.perm), w.perm^(-1))
end

