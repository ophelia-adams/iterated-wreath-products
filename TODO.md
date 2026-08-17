# TODO

## permutations.jl
- 

## wreath.jl
- Make `IteratedWreath{N}`'s equality work (`nothing` should be the same as the identity in `Perm{N}`).

## wreathexpr.jl
- Figure out how to match on parametrized types in the rules to eliminate helper functions `wreathdim` and `wreathtype`.
- How to make rules associative etc. (is flattening the way?).
- Maybe adding exponent type to group things like `x * x * x`?
- `Base.show` for these awful expression types.

## Other/general
- Macros or something to make all this more usable. (low priority)
- Substitution dictionaries for mutual recursions.
- Computations on elements.
