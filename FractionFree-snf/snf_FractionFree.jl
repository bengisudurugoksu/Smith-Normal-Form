module FractionFreeSNF_Unimodular

using LinearAlgebra

export extended_gcd, snf_fraction_free_unimodular, verify_snf

# --------------------------------------------------
# Extended Euclidean Algorithm: returns (g,x,y) with x*a + y*b = g >= 0
# --------------------------------------------------
function extended_gcd(a::T, b::T) where {T<:Integer}
    if b == 0
        return (abs(a), sign(a), zero(T))
    end
    old_r, r = a, b
    old_s, s = one(T), zero(T)
    old_t, t = zero(T), one(T)
    while r != 0
        q = div(old_r, r)
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    end
    old_r < 0 ? (-old_r, -old_s, -old_t) : (old_r, old_s, old_t)
end
extended_gcd(a::Integer, b::Integer) = extended_gcd(promote(a,b)...)

# --------------------------------------------------
# Helpers
# --------------------------------------------------
@inline function swap_rows!(A, i::Int, j::Int)
    i == j && return
    A[i, :], A[j, :] = A[j, :], copy(A[i, :])
end

@inline function swap_cols!(A, i::Int, j::Int)
    i == j && return
    A[:, i], A[:, j] = A[:, j], copy(A[:, i])
end

# Find pivot with smallest abs value in submatrix (k:m, k:n)
function find_best_pivot(S, k::Int)
    m, n = size(S)
    best_i, best_j = 0, 0
    best_abs = typemax(eltype(S))
    for i in k:m
        for j in k:n
            v = S[i,j]
            if v != 0
                av = abs(v)
                if av < best_abs
                    best_abs = av
                    best_i, best_j = i, j
                end
            end
        end
    end
    return best_i, best_j
end

# --------------------------------------------------
# Unimodular gcd-based elimination on rows:
# For a = S[k,k], b = S[i,k], produce new (a,b) -> (g,0) in pivot column
# and update U accordingly so U*A*V = S always holds.
# --------------------------------------------------
function elim_row_unimod!(S, U, k::Int, i::Int)
    b = S[i,k]
    b == 0 && return false
    a = S[k,k]
    g, x, y = extended_gcd(a, b)
    u = div(a, g)  # exact
    v = div(b, g)  # exact

    Sk = copy(S[k,:]); Si = copy(S[i,:])
    Uk = copy(U[k,:]); Ui = copy(U[i,:])

    S[k,:] =  x .* Sk .+ y .* Si
    S[i,:] = -v .* Sk .+ u .* Si

    U[k,:] =  x .* Uk .+ y .* Ui
    U[i,:] = -v .* Uk .+ u .* Ui
    return true
end

# --------------------------------------------------
# Unimodular gcd-based elimination on columns:
# For a = S[k,k], c = S[k,j], produce new (a,c) -> (g,0) in pivot row
# and update V accordingly.
# --------------------------------------------------
function elim_col_unimod!(S, V, k::Int, j::Int)
    c = S[k,j]
    c == 0 && return false
    a = S[k,k]
    g, x, y = extended_gcd(a, c)
    u = div(a, g)
    v = div(c, g)

    Ck = copy(S[:,k]); Cj = copy(S[:,j])
    Vk = copy(V[:,k]); Vj = copy(V[:,j])

    S[:,k] =  x .* Ck .+ y .* Cj
    S[:,j] = -v .* Ck .+ u .* Cj

    V[:,k] =  x .* Vk .+ y .* Vj
    V[:,j] = -v .* Vk .+ u .* Vj
    return true
end

# pivot divides everything in lower-right block?
function find_div_violation(S, k::Int)
    m, n = size(S)
    d = S[k,k]
    d == 0 && return nothing
    for i in k+1:m
        for j in k+1:n
            if S[i,j] != 0 && rem(S[i,j], d) != 0
                return (i,j)
            end
        end
    end
    return nothing
end

# Inject violating entry into pivot column/row using unimodular add:
# col_k += col_vj (or row_k += row_vi) keeps V/U unimodular.
@inline function add_col!(S, V, k::Int, j::Int)
    S[:,k] .+= S[:,j]
    V[:,k] .+= V[:,j]
end

@inline function add_row!(S, U, k::Int, i::Int)
    S[k,:] .+= S[i,:]
    U[k,:] .+= U[i,:]
end

# --------------------------------------------------
# Main SNF (fraction-free in the gcd-exact sense), returns S,U,V unimodular.
# --------------------------------------------------
function snf_fraction_free_unimodular(A::AbstractMatrix{T}) where {T<:Integer}
    m, n = size(A)
    S = copy(A)
    U = Matrix{T}(I, m, m)
    V = Matrix{T}(I, n, n)
    r = min(m,n)

    # -------- Phase 1: Diagonalize with "pivot divides remaining block" invariant --------
    for k in 1:r
        pi, pj = find_best_pivot(S, k)
        pi == 0 && break

        swap_rows!(S, k, pi); swap_rows!(U, k, pi)
        swap_cols!(S, k, pj); swap_cols!(V, k, pj)

        # Iterate until cleared and divisibility-in-block holds
        guard = 0
        while true
            guard += 1
            guard > 8000 && error("Guard triggered: potential non-convergence at k=$k")

            changed = false

            # Clear pivot column fully (all i != k)
            for i in 1:m
                i == k && continue
                changed |= elim_row_unimod!(S, U, k, i)
            end

            # Clear pivot row fully (all j != k)
            for j in 1:n
                j == k && continue
                changed |= elim_col_unimod!(S, V, k, j)
            end

            # Make pivot positive
            if S[k,k] < 0
                S[k,:] .*= -1
                U[k,:] .*= -1
            end

            viol = find_div_violation(S, k)
            if viol === nothing
                break
            else
                vi, vj = viol
                # Unimodular "injection" to force gcd reduction next iteration
                add_col!(S, V, k, vj)
                changed = true
            end

            !changed && break
        end
    end

    # -------- Phase 2: Normalize invariant factors (enforce d_k | d_{k+1}) --------
    # We repair adjacent diagonal pairs; after each repair we re-clear to restore diagonality.
    for k in 1:(r-1)
        guard = 0
        while S[k,k] != 0 && S[k+1,k+1] != 0 && rem(S[k+1,k+1], S[k,k]) != 0
            guard += 1
            guard > 2000 && error("Normalization guard triggered at k=$k")

            a = S[k,k]
            b = S[k+1,k+1]
            g, x, y = extended_gcd(a, b)
            u = div(a, g)
            v = div(b, g)

            # Apply unimodular transform on rows k,k+1
            Rk = copy(S[k,:]); R2 = copy(S[k+1,:])
            Uk = copy(U[k,:]); U2 = copy(U[k+1,:])
            S[k,:]   =  x .* Rk .+ y .* R2
            S[k+1,:] = -v .* Rk .+ u .* R2
            U[k,:]   =  x .* Uk .+ y .* U2
            U[k+1,:] = -v .* Uk .+ u .* U2

            # Apply unimodular transform on cols k,k+1
            Ck = copy(S[:,k]); C2 = copy(S[:,k+1])
            Vk = copy(V[:,k]); V2 = copy(V[:,k+1])
            S[:,k]   =  x .* Ck .+ y .* C2
            S[:,k+1] = -v .* Ck .+ u .* C2
            V[:,k]   =  x .* Vk .+ y .* V2
            V[:,k+1] = -v .* Vk .+ u .* V2

            # Re-clear around pivots k and k+1 to restore diagonality
            for i in 1:m
                i == k && continue
                elim_row_unimod!(S, U, k, i)
            end
            for i in 1:m
                i == k+1 && continue
                elim_row_unimod!(S, U, k+1, i)
            end
            for j in 1:n
                j == k && continue
                elim_col_unimod!(S, V, k, j)
            end
            for j in 1:n
                j == k+1 && continue
                elim_col_unimod!(S, V, k+1, j)
            end

            if S[k,k] < 0
                S[k,:] .*= -1
                U[k,:] .*= -1
            end
        end
    end

    # final sign normalization
    for k in 1:r
        if S[k,k] < 0
            S[k,:] .*= -1
            U[k,:] .*= -1
        end
    end

    return S, U, V
end

# --------------------------------------------------
# Verification (recommended)
# --------------------------------------------------
function verify_snf(A, S, U, V)
    if U*A*V != S
        return false
    end
    # diagonal check
    m, n = size(S)
    for i in 1:m, j in 1:n
        if i != j && S[i,j] != 0
            return false
        end
    end
    # divisibility check
    d = [S[i,i] for i in 1:min(m,n) if S[i,i] != 0]
    for i in 1:(length(d)-1)
        if rem(d[i+1], d[i]) != 0
            return false
        end
    end
    return true
end

end # module
