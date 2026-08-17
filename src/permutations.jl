# Original code: https://github.com/ophelia-adams/Monodromy.jl/blob/main/src/core/permutations.jl
# Original Author: Ophelia Adams
# ================================

"""
An abstract type that encompasses both `Σ{N}` and `Cyc{N}`. Used in wreath products.
"""
abstract type Perm{N} end


# ======== Permutation Groups ======== #


"""
    Σ{N}

Wrapper over an N-tuple representing a permutation. The inner constructors check it is a permutation.
"""
struct Σ{N} <: Perm{N}
    perm::NTuple{N, Int}
    name::Union{String, Nothing}

    function Σ{N}(p::NTuple{N, Int}, name=nothing) where N
        Set(p) == Set(1:N) || throw("$p is not a valid permutation.")
        new{N}(p, name)
    end

    function Σ(p::NTuple{N, Int}, name=nothing) where N
        Set(p) == Set(1:N) || throw("$p is not a valid permutation.")
        new{N}(p, name)
    end

    # Variant of constructor that returns the identity element.
    Σ{N}() where N = new{N}(Tuple(1:N), nothing)
end

# identity permutation
Base.one(::Σ{N}) where N = Σ{N}()
Base.one(::Type{Σ{N}}) where N = Σ{N}()

Base.getindex(p::Σ{N}, i::Int) where N = p.perm[i]
Base.isempty(p::Σ{N}) where N = isempty(p.perm)
# Base.keys(p::Σ{N})::Vector{Symbol} = [k for k in keys(p.perm)]
# Base.values(p::Σ{N})::Vector{Symbol} = [v for v in values(p.perm)]

function Base.show(io::IO, p::Σ{N}) where N
    if p.name == nothing
        print(io, cyclestring(p))
    else
        print(io, p.name)
    end
end

"""
    cycleof(p::Σ{N}, i::Int)::Vector{Int}

Calculates the cycle, or orbit, of permutation `p` on symbol `i`.
"""
function cycleof(p::Σ{N}, i::Int)::Vector{Int} where N
    cycle = Vector{Int}()
    push!(cycle, i)
    while !(p[cycle[end]] == i)
        push!(cycle, p[cycle[end]])
    end
    return cycle
end

"""
    cycles(p::Σ{N})::Vector{Vector{Int}}

Returns a representation of `p` as a list of disjoint cycles.
"""
function cycles(p::Σ{N}) :: Vector{Vector{Int}} where N
    remaining = [N-i+1 for i in 1:N]
    cycles = Vector{Vector{Int}}()

    while !(isempty(remaining))
        name = pop!(remaining)
        cycle = cycleof(p, name)
        push!(cycles, cycle)
        filter!(x -> x ∉ cycle, remaining)
    end

    return cycles
end

"""
    cyclestring(p::Σ{N})

Returns a string with the decomposition of `p` into disjoint cycles.
"""
function cyclestring(p::Σ{N}) where N
    cycs = cycles(p)
    s = ""
    for c in cycs
        if length(c) > 1
            s = s * "("
            s = s * prod([string(l) * ", " for l in c[1:end-1]])
            s = s * string(c[end]) * ")"
        end
    end
    s = s == "" ? "1" : s
    return s
end

"""
    fromcycle(cycle::Vector{Int}, N::Int)::Σ{N}

Returns a permutation given a valid cycle.
"""
function fromcycle(cycle::Vector{Int}, N::Int)::Σ{N}
    if isempty(cycle)
        return Σ{N}()
    end
    valid = length(cycle) == length(unique(cycle)) <= N
    valid = valid && maximum(cycle) <= N
    valid || throw("$cycle is not a valid cycle in Σ{$N}.")
    cyc = [cycle; cycle[1]]
    p = Tuple(if (k = findfirst(==(i), cyc)) !== nothing
                  cyc[k+1]
              else
                  i
              end for i in 1:N)
    Σ{N}(p)
end

"""
    fromcycles(cycles::Vector{Vector{Int}})::Σ{N}

Returns the permutation specified by the list of disjoint cycles.
"""
function fromcycles(cycles::Vector{Vector{Int}}, N::Int)::Σ{N}
    cycles1 = filter(!isempty, cycles)
    catcycles = vcat(cycles1...)
    if isempty(catcycles)
        return Σ{N}()
    end
    valid = length(catcycles) == length(unique(catcycles)) <= N
    valid = valid && maximum(catcycles) <= N
    valid || throw("$cycles is not a valid collection of disjoint cycles in Σ{$N}.")
    cycs = map(x -> [x; x[1]], cycles1)
    catcycs = vcat(cycs...)
    p = Tuple(if (k = findfirst(==(i), catcycs)) !== nothing
                  catcycs[k+1]
              else
                  i
              end for i in 1:N)
    Σ{N}(p)
end

"""
    inverse(p::Σ{N})::Σ{N}

Returns the inverse of a permutation.
"""
function inverse(p::Σ{N})::Σ{N} where N
    cycs = cycles(p)
    revcycs = map(reverse, cycs)
    fromcycles(revcycs, N)
end

Base.inv(p::Σ{N}) where N = inverse(p)

"""
    ((s::Int)^(p::Σ{N}))::Int

Applies the permutation `p` to `s`, acting on the *right*.
"""
Base.:^(s::Int, p::Σ{N}) where N = p[s]

"""
    ((p::Σ{N})*(q::Σ{N}))::Σ{N}

Composes permutations `p` and `q`, using the more natural *right* action, so that
    a^(p*q) = (a^p)^q
"""
function Base.:*(p::Σ{N}, q::Σ{N})::Σ{N} where N
    Σ{N}(Tuple(q[p[i]] for i in 1:N))
end

"""
    ((p::Σ{N})^(q::Σ{N}))::Σ{N}

Conjugates the permutation `p` by the permutation `q`.
"""
Base.:^(p::Σ{N}, q::Σ{N}) where N = inverse(q)*p*q

"""
    ((p::Σ{N})^(n::Int))::Σ{N}

Raises `p` to the `n`-th power.
"""
function Base.:^(p::Σ{N}, n::Int)::Σ{N} where N
    if n == 0
        Σ{N}()
    elseif n > 0
        cycs = cycles(p)
        advance_n = c ->
            if mod(n, length(c)) == 0
                []
            else
                [c[mod1(1 + (i-1)*n, length(c))] for i in 1:length(c)]
            end
        cycs_shifted = map(advance_n, cycs)
        if isempty(vcat(cycs_shifted...))
            Σ{N}()
        else
            fromcycles(cycs_shifted, N)
        end
    else
        inverse(p)^(-n)
    end
end


# ======== Cyclic Groups ======== #


"""
    Cyc{N}

Isomorphic to `Z/NZ`. Supposed to represent the subgroup generated by the `N`-cycle `(1 2 ... N)`.
"""
struct Cyc{N} <: Perm{N}
    pow::Int
    name::String

    function Cyc{N}(k::Int, name="σ") :: Cyc{N} where N
        new{N}(mod(k, N), name)
    end

    Cyc{N}() where N = Cyc{N}(0)
end

# identity element
Base.one(::Cyc{N}) where N = Cyc{N}()
Base.one(::Type{Cyc{N}}) where N = Cyc{N}()

Base.getindex(g::Cyc{N}, i::Int) where N = mod1(i + g.pow, N)

function Base.show(io::IO, g::Cyc{N}) where N
    if g.pow == 0
        print(io, "1")
    elseif g.pow == 1
        print(io, g.name)
    else
        print(io, g.name * superscript(g.pow))
    end
end

# Yes, this is awful (and ChatGPT'd). Alas, it works.
superscript_digits = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']

# Helper function for displaying exponents.
function superscript(k::Integer)
    k < 0 && return "⁻" * superscript(-k)
    join(superscript_digits[d + 1] for d in reverse(digits(k, base=10)))
end

Base.inv(g::Cyc{N}) where N = Cyc{N}(-g.pow, g.name)

"""
    ((s::Int)^(g::Cyc{N}))::Int

Applies `g` to `s` as a cyclic permutation.
"""
Base.:^(s::Int, g::Cyc{N}) where N = g[s]

"""
    ((g::Cyc{N})*(h::Cyc{N}))::Σ{N}

Multiplies `g` and `h` in the cyclic group of `N` elements.
Inherits its string representation from the first argument.
"""
function Base.:*(g::Cyc{N}, h::Cyc{N})::Cyc{N} where N
    Cyc{N}(g.pow + h.pow, g.name)
end

"""
    ((g::Cyc{N})^(n::Int))::Cyc{N}

Raises `g` to the `n`-th power. 
"""
function Base.:^(g::Cyc{N}, n::Int)::Cyc{N} where N
    Cyc{N}(g.pow * n, g.name)
end

