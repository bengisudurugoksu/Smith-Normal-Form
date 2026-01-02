module ModularNemoSNF

using LinearAlgebra
using Nemo

export snf_modular_nemo

const ZZ = Nemo.ZZ

"""
    snf_modular_nemo(A::AbstractMatrix{<:Integer})

Computes Smith Normal Form using Nemo.jl.
Internally uses modular algorithms + CRT.
Returns (S, U, V) where U and V are identity matrices.
"""
function snf_modular_nemo(A::AbstractMatrix{<:Integer})
    m, n = size(A)

    A_ZZ = matrix(ZZ, m, n, vec(A))
    S_ZZ = snf(A_ZZ)

    S = Matrix{Int}(undef, m, n)
    for i in 1:m, j in 1:n
        S[i,j] = Int(S_ZZ[i,j])
    end

    U = Matrix{Int}(I, m, m)
    V = Matrix{Int}(I, n, n)

    return S, U, V
end

end # module
