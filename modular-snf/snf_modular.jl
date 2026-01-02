# Modular Smith Normal Form (Diagonal only, no U/V)
# - No external packages
# - Uses modular determinants + CRT reconstruction
# - Computes determinantal divisors Δ_k = gcd of all k×k minors
# - Invariant factors: s1=Δ1, sk=Δk/Δ(k-1)
############################################################

module ModularSNF

using LinearAlgebra

export modular_snf

# ----------------------------
# Basic number theory helpers
# ----------------------------

# Extended GCD: returns (g, x, y) such that a*x + b*y = g = gcd(a,b)
function egcd(a::Int, b::Int)
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0
        q = div(old_r, r)
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    end
    # Normalize gcd to be positive
    if old_r < 0
        return (-old_r, -old_s, -old_t)
    end
    return (old_r, old_s, old_t)
end

# Modular inverse of a mod p (p prime, but works if gcd(a,p)=1)
function invmod(a::Int, p::Int)
    a = mod(a, p)
    g, x, _ = egcd(a, p)
    if g != 1
        error("invmod: inverse does not exist (gcd != 1).")
    end
    return mod(x, p)
end

# Combine x ≡ a (mod m) and x ≡ b (mod p) into x ≡ a2 (mod m*p)
# Assumes gcd(m,p)=1
function crt_pair(a::Int, m::Int, b::Int, p::Int)
    # Solve: a + m*t ≡ b (mod p)  =>  m*t ≡ (b-a) (mod p)
    t = mod((b - a), p)
    inv_m = invmod(mod(m, p), p)
    t = mod(t * inv_m, p)
    a2 = a + m * t
    m2 = m * p
    a2 = mod(a2, m2)
    return a2, m2
end

# Symmetric representative in (-M/2, M/2]
function symm_rep(x::Int, M::Int)
    x = mod(x, M)
    if x > div(M, 2)
        x -= M
    end
    return x
end

# gcd for Int with absolute safety
gcd_abs(a::Int, b::Int) = gcd(abs(a), abs(b))

# ----------------------------
# Prime generator (simple)
# ----------------------------
function isprime(n::Int)
    if n < 2
        return false
    end
    if n % 2 == 0
        return n == 2
    end
    d = 3
    while d*d <= n
        if n % d == 0
            return false
        end
        d += 2
    end
    return true
end

function first_primes(count::Int; start::Int=3)
    ps = Int[]
    x = start
    while length(ps) < count
        if isprime(x)
            push!(ps, x)
        end
        x += 1
    end
    return ps
end

# ----------------------------
# Determinant mod p
# Gaussian elimination over Z/pZ
# ----------------------------
function det_mod_p(A::Matrix{Int}, p::Int)
    n = size(A, 1)
    @assert size(A, 2) == n

    M = Array{Int}(undef, n, n)
    for i in 1:n, j in 1:n
        M[i,j] = mod(A[i,j], p)
    end

    det = 1
    sign = 1

    for col in 1:n
        pivot = 0
        pivot_row = 0
        for r in col:n
            if M[r,col] != 0
                pivot = M[r,col]
                pivot_row = r
                break
            end
        end
        if pivot == 0
            return 0
        end

        if pivot_row != col
            # swap rows
            for j in col:n
                M[col,j], M[pivot_row,j] = M[pivot_row,j], M[col,j]
            end
            sign = -sign
            pivot = M[col,col]
        end

        det = mod(det * pivot, p)
        inv_pivot = invmod(pivot, p)

        # eliminate rows below
        for r in (col+1):n
            if M[r,col] != 0
                factor = mod(M[r,col] * inv_pivot, p)
                for j in col:n
                    M[r,j] = mod(M[r,j] - factor * M[col,j], p)
                end
            end
        end
    end

    det = mod(det * (sign == -1 ? (p-1) : 1), p)
    return det
end

# ----------------------------
# Hadamard bound for determinant of an integer matrix
# |det(A)| <= Π ||row_i||_2
# This is safe for deciding when CRT modulus is "large enough".
# ----------------------------
function hadamard_bound(A::Matrix{Int})
    n = size(A, 1)
    @assert size(A, 2) == n
    prod_bound = 1.0
    for i in 1:n
        s = 0.0
        for j in 1:n
            s += float(A[i,j])^2
        end
        prod_bound *= sqrt(s)
    end
    # ceil to Int; also handle near-zero
    b = Int(ceil(prod_bound))
    return max(b, 1)
end

# ----------------------------
# Extract submatrix given row/col index sets
# ----------------------------
function submatrix(A::Matrix{Int}, rows::Vector{Int}, cols::Vector{Int})
    k = length(rows)
    B = Array{Int}(undef, k, k)
    for i in 1:k, j in 1:k
        B[i,j] = A[rows[i], cols[j]]
    end
    return B
end

# ----------------------------
# Generate combinations of indices (1..n choose k)
# ----------------------------
function combinations(n::Int, k::Int)
    result = Vector{Vector{Int}}()
    comb = collect(1:k)

    function push_current()
        push!(result, copy(comb))
    end

    if k == 0
        push!(result, Int[])
        return result
    end
    if k > n
        return result
    end

    push_current()
    while true
        i = k
        while i >= 1 && comb[i] == n - k + i
            i -= 1
        end
        if i == 0
            break
        end
        comb[i] += 1
        for j in (i+1):k
            comb[j] = comb[j-1] + 1
        end
        push_current()
    end
    return result
end

# ----------------------------
# Reconstruct exact determinant using modular primes + CRT
# stops when modulus M > 2*HadamardBound
# ----------------------------
function det_integer_via_crt(A::Matrix{Int}; primes::Vector{Int}=first_primes(50))
    n = size(A, 1)
    @assert size(A, 2) == n

    bound = hadamard_bound(A)
    target = 2 * bound + 1

    a = 0      # current residue
    M = 1      # current modulus

    for p in primes
        d = det_mod_p(A, p)
        a, M = crt_pair(a, M, d, p)
        if M >= target
            return symm_rep(a, M)
        end
    end

    # If we run out of primes, still return best symmetric rep we have
    return symm_rep(a, M)
end

# ----------------------------
# Compute Δ_k: gcd of all k×k minors' determinants
# (determinants reconstructed via CRT)
# ----------------------------
function determinantal_divisor(A::Matrix{Int}, k::Int; primes::Vector{Int}=first_primes(80))
    m, n = size(A)
    @assert 1 <= k <= min(m,n)

    row_combs = combinations(m, k)
    col_combs = combinations(n, k)

    g = 0
    for rows in row_combs
        for cols in col_combs
            B = submatrix(A, rows, cols)
            d = det_integer_via_crt(B; primes=primes)
            g = (g == 0) ? abs(d) : gcd_abs(g, d)
            # Early exit: gcd becomes 1, cannot get smaller
            if g == 1
                return 1
            end
        end
    end
    return g
end

# ----------------------------
# Modular SNF diagonal (invariant factors only)
# Returns a vector diag of length min(m,n) with invariant factors
# Remaining (if m!=n) can be interpreted as zeros beyond rank.
# ----------------------------
function modular_snf_diagonal(A::Matrix{Int}; primes::Vector{Int}=first_primes(120))
    m, n = size(A)
    rmax = min(m, n)

    Δ_prev = 1
    diag = Int[]
    rank = 0

    for k in 1:rmax
        Δk = determinantal_divisor(A, k; primes=primes)
        if Δk == 0
            break
        end
        # Invariant factor s_k = Δ_k / Δ_{k-1}
        sk = div(Δk, Δ_prev)
        push!(diag, sk)
        Δ_prev = Δk
        rank = k
    end

    # Fill remaining with zeros (optional; helpful for "diagonal SNF matrix")
    while length(diag) < rmax
        push!(diag, 0)
    end

    return diag
end

# ----------------------------
# Build the diagonal SNF matrix from diag vector
# ----------------------------
function diag_matrix(diag::Vector{Int}, m::Int, n::Int)
    D = fill(0, m, n)
    for i in 1:min(length(diag), min(m,n))
        D[i,i] = diag[i]
    end
    return D
end

############################################################
# Example usage
############################################################
# A = [ 2 4 4;
#       6 6 12;
#       10 4 16 ]
#
# diag = modular_snf_diagonal(A)
# D = diag_matrix(diag, size(A,1), size(A,2))
# println("SNF diagonal invariants: ", diag)
# println("SNF diagonal matrix:\n", D)
############################################################

# ----------------------------
# Wrapper for compatibility with other algorithms
# Returns (S, U, V) where S is diagonal SNF matrix
# U and V are identity matrices (not computed by modular method)
# ----------------------------
function modular_snf(A::AbstractMatrix{T}) where T <: Integer
    m, n = size(A)
    A_int = Matrix{Int}(A)
    diag = modular_snf_diagonal(A_int)
    S = diag_matrix(diag, m, n)
    U = Matrix{T}(I, m, m)
    V = Matrix{T}(I, n, n)
    return S, U, V
end

end # module ModularSNF
