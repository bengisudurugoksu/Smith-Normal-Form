"""
Smith Normal Form Implementation using Extended Euclidean Algorithm
===================================================================

This module computes the Smith Normal Form (SNF) for integer matrices.
Unimodular transformations are obtained using the Extended Euclidean Algorithm (EEA).

Theory:
------
For every m×n integer matrix A, there exist unimodular matrices U (m×m) and V (n×n) 
such that:
    S = U * A * V
where S is a diagonal matrix and its diagonal elements d₁, d₂, ..., dᵣ
(r = rank(A)) satisfy: d₁ | d₂ | ... | dᵣ (each element divides the next)

References:
-----------
1. Cohen, H. (1993). A Course in Computational Algebraic Number Theory. Springer.
2. Storjohann, A. (2000). Algorithms for Matrix Canonical Forms. PhD Thesis.
"""

module SmithNormalForm

using LinearAlgebra

export extended_gcd, smith_normal_form, verify_snf, elementary_divisors

"""
    extended_gcd(a::Integer, b::Integer) -> (gcd, x, y)

Extended Euclidean Algorithm.
Finds Bézout coefficients: a*x + b*y = gcd(a, b)

# Example
```julia
julia> g, x, y = extended_gcd(35, 15)
(5, 1, -2)
julia> 35*x + 15*y == g
true
```
"""
function extended_gcd(a::T, b::T) where T <: Integer
    if b == 0
        return (abs(a), sign(a), zero(T))
    end
    
    old_r, r = a, b
    old_s, s = one(T), zero(T)
    old_t, t = zero(T), one(T)
    
    while r != 0
        quotient = div(old_r, r)
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    end
    
    # GCD must be positive
    if old_r < 0
        return (-old_r, -old_s, -old_t)
    end
    
    return (old_r, old_s, old_t)
end

# For different type arguments
extended_gcd(a::Integer, b::Integer) = extended_gcd(promote(a, b)...)

"""
    swap_rows!(A, i, j)

Swaps rows i and j.
"""
function swap_rows!(A::AbstractMatrix, i::Int, j::Int)
    if i != j
        A[i, :], A[j, :] = A[j, :], copy(A[i, :])
    end
end

"""
    swap_cols!(A, i, j)

Swaps columns i and j.
"""
function swap_cols!(A::AbstractMatrix, i::Int, j::Int)
    if i != j
        A[:, i], A[:, j] = A[:, j], copy(A[:, i])
    end
end

"""
    add_row_multiple!(A, target, source, factor)

A[target, :] += factor * A[source, :]
"""
function add_row_multiple!(A::AbstractMatrix, target::Int, source::Int, factor)
    A[target, :] .+= factor .* A[source, :]
end

"""
    add_col_multiple!(A, target, source, factor)

A[:, target] += factor * A[:, source]
"""
function add_col_multiple!(A::AbstractMatrix, target::Int, source::Int, factor)
    A[:, target] .+= factor .* A[:, source]
end

"""
    find_pivot(A, start_row, start_col)

Finds the smallest non-zero element (in absolute value) 
starting from start_row and start_col.
"""
function find_pivot(A::AbstractMatrix, start_row::Int, start_col::Int)
    m, n = size(A)
    min_val = typemax(eltype(A))
    min_i, min_j = 0, 0
    
    for i in start_row:m
        for j in start_col:n
            if A[i, j] != 0 && abs(A[i, j]) < min_val
                min_val = abs(A[i, j])
                min_i, min_j = i, j
            end
        end
    end
    
    return min_i, min_j
end

"""
    eliminate_row!(S, U, pivot_row, pivot_col)

Eliminates other rows in the pivot column using EEA.
"""
function eliminate_row!(S::AbstractMatrix, U::AbstractMatrix, 
                        pivot_row::Int, pivot_col::Int)
    m = size(S, 1)
    changed = false
    
    for i in 1:m
        if i != pivot_row && S[i, pivot_col] != 0
            a = S[pivot_row, pivot_col]
            b = S[i, pivot_col]
            g, x, y = extended_gcd(a, b)
            
            # Unimodular transformation
            u = div(a, g)
            v = div(b, g)
            
            # New rows
            new_pivot_row = x .* S[pivot_row, :] .+ y .* S[i, :]
            new_i_row = (-v) .* S[pivot_row, :] .+ u .* S[i, :]
            
            S[pivot_row, :] = new_pivot_row
            S[i, :] = new_i_row
            
            # Update U matrix
            new_pivot_U = x .* U[pivot_row, :] .+ y .* U[i, :]
            new_i_U = (-v) .* U[pivot_row, :] .+ u .* U[i, :]
            
            U[pivot_row, :] = new_pivot_U
            U[i, :] = new_i_U
            
            changed = true
        end
    end
    
    return changed
end

"""
    eliminate_col!(S, V, pivot_row, pivot_col)

Eliminates other columns in the pivot row using EEA.
"""
function eliminate_col!(S::AbstractMatrix, V::AbstractMatrix,
                        pivot_row::Int, pivot_col::Int)
    n = size(S, 2)
    changed = false
    
    for j in 1:n
        if j != pivot_col && S[pivot_row, j] != 0
            a = S[pivot_row, pivot_col]
            b = S[pivot_row, j]
            g, x, y = extended_gcd(a, b)
            
            u = div(a, g)
            v = div(b, g)
            
            # New columns
            new_pivot_col = x .* S[:, pivot_col] .+ y .* S[:, j]
            new_j_col = (-v) .* S[:, pivot_col] .+ u .* S[:, j]
            
            S[:, pivot_col] = new_pivot_col
            S[:, j] = new_j_col
            
            # Update V matrix
            new_pivot_V = x .* V[:, pivot_col] .+ y .* V[:, j]
            new_j_V = (-v) .* V[:, pivot_col] .+ u .* V[:, j]
            
            V[:, pivot_col] = new_pivot_V
            V[:, j] = new_j_V
            
            changed = true
        end
    end
    
    return changed
end

"""
    check_divisibility(S, pivot_row, pivot_col)

Checks whether the pivot element divides all elements 
in the lower-right submatrix.
"""
function check_divisibility(S::AbstractMatrix, pivot_row::Int, pivot_col::Int)
    m, n = size(S)
    d = S[pivot_row, pivot_col]
    
    if d == 0
        return true, 0, 0
    end
    
    for i in (pivot_row + 1):m
        for j in (pivot_col + 1):n
            if S[i, j] != 0 && S[i, j] % d != 0
                return false, i, j
            end
        end
    end
    
    return true, 0, 0
end

"""
    smith_normal_form(A::AbstractMatrix{<:Integer})

Computes the Smith Normal Form of integer matrix A.

# Returns
- `S`: Smith Normal Form (diagonal matrix)
- `U`: Left unimodular matrix (det(U) = ±1)
- `V`: Right unimodular matrix (det(V) = ±1)

# Relation
U * A * V = S

# Example
```julia
julia> A = [2 4 4; -6 6 12; 10 4 16]
julia> S, U, V = smith_normal_form(A)
julia> U * A * V == S
true
```
"""
function smith_normal_form(A::AbstractMatrix{T}) where T <: Integer
    m, n = size(A)
    
    # Working copies
    S = copy(A)
    U = Matrix{T}(I, m, m)  # m×m identity matrix
    V = Matrix{T}(I, n, n)  # n×n identity matrix
    
    min_dim = min(m, n)
    
    for k in 1:min_dim
        # Find pivot
        pi, pj = find_pivot(S, k, k)
        
        if pi == 0  # Remaining submatrix is all zeros
            break
        end
        
        # Move pivot to position (k, k)
        if pi != k
            swap_rows!(S, k, pi)
            swap_rows!(U, k, pi)
        end
        if pj != k
            swap_cols!(S, k, pj)
            swap_cols!(V, k, pj)
        end
        
        # Row and column elimination to zero out elements except (k, k)
        # Loop until divisibility is satisfied
        max_iterations = 1000  # Infinite loop protection
        iteration = 0
        
        while iteration < max_iterations
            iteration += 1
            
            # First eliminate
            while true
                row_changed = eliminate_row!(S, U, k, k)
                col_changed = eliminate_col!(S, V, k, k)
                
                if !row_changed && !col_changed
                    break
                end
            end
            
            # Divisibility check - for entire lower-right block
            divisible, bad_i, bad_j = check_divisibility(S, k, k)
            
            if divisible
                break  # OK, move to next k
            else
                # Add bad_j column to pivot column, then eliminate again
                add_col_multiple!(S, k, bad_j, one(T))
                add_col_multiple!(V, k, bad_j, one(T))
            end
        end
        
        # Make diagonal element positive
        if S[k, k] < 0
            S[k, :] .*= -1
            U[k, :] .*= -1
        end
    end
    
    # Final divisibility correction - ensure d[i] | d[i+1]
    for i in 1:(min_dim - 1)
        if S[i, i] != 0 && S[i+1, i+1] != 0 && S[i+1, i+1] % S[i, i] != 0
            # Find GCD of d[i] and d[i+1] and correct
            g, x, y = extended_gcd(S[i, i], S[i+1, i+1])
            u = div(S[i, i], g)
            v = div(S[i+1, i+1], g)
            
            # Row transformation
            new_row_i = x .* S[i, :] .+ y .* S[i+1, :]
            new_row_i1 = (-v) .* S[i, :] .+ u .* S[i+1, :]
            S[i, :] = new_row_i
            S[i+1, :] = new_row_i1
            
            new_U_i = x .* U[i, :] .+ y .* U[i+1, :]
            new_U_i1 = (-v) .* U[i, :] .+ u .* U[i+1, :]
            U[i, :] = new_U_i
            U[i+1, :] = new_U_i1
            
            # Column transformation to make diagonal
            if S[i, i+1] != 0
                factor = -div(S[i, i+1], S[i, i])
                add_col_multiple!(S, i+1, i, factor)
                add_col_multiple!(V, i+1, i, factor)
            end
            if S[i+1, i] != 0
                factor = -div(S[i+1, i], S[i+1, i+1])
                add_row_multiple!(S, i, i+1, factor)
                add_row_multiple!(U, i, i+1, factor)
            end
        end
    end
    
    # Make positive
    for k in 1:min_dim
        if S[k, k] < 0
            S[k, :] .*= -1
            U[k, :] .*= -1
        end
    end
    
    return S, U, V
end

"""
    verify_snf(A, S, U, V) -> Bool

Verifies the correctness of SNF computation.

Checks:
1. U * A * V = S
2. det(U) = ±1 (unimodular)
3. det(V) = ±1 (unimodular)
4. S is a diagonal matrix
5. d[i] | d[i+1] (divisibility)
"""
function verify_snf(A::AbstractMatrix, S::AbstractMatrix, 
                    U::AbstractMatrix, V::AbstractMatrix)
    m, n = size(A)
    
    # 1. U * A * V = S check
    if U * A * V != S
        println("❌ U * A * V ≠ S")
        return false
    end
    println("✓ U * A * V = S")
    
    # 2. Unimodular check (det = ±1)
    det_U = Int(round(det(U)))
    det_V = Int(round(det(V)))
    
    if abs(det_U) != 1
        println("❌ det(U) = $det_U ≠ ±1")
        return false
    end
    println("✓ det(U) = $det_U")
    
    if abs(det_V) != 1
        println("❌ det(V) = $det_V ≠ ±1")
        return false
    end
    println("✓ det(V) = $det_V")
    
    # 3. S diagonal matrix check
    for i in 1:m
        for j in 1:n
            if i != j && S[i, j] != 0
                println("❌ S is not diagonal: S[$i,$j] = $(S[i,j])")
                return false
            end
        end
    end
    println("✓ S is diagonal")
    
    # 4. Divisibility check (d[i] | d[i+1])
    diag = [S[i,i] for i in 1:min(m,n) if S[i,i] != 0]
    for i in 1:(length(diag)-1)
        if diag[i+1] % diag[i] != 0
            println("❌ d[$i] = $(diag[i]) does not divide d[$(i+1)] = $(diag[i+1])")
            return false
        end
    end
    println("✓ Divisibility condition satisfied")
    
    println("\n✅ All verification tests passed!")
    return true
end

"""
    elementary_divisors(A::AbstractMatrix{<:Integer})

Returns the elementary divisors (invariant factors) of the matrix.
These are the non-zero diagonal elements of the SNF.
"""
function elementary_divisors(A::AbstractMatrix{T}) where T <: Integer
    S, _, _ = smith_normal_form(A)
    m, n = size(S)
    return [S[i,i] for i in 1:min(m,n) if S[i,i] != 0]
end

end # module
