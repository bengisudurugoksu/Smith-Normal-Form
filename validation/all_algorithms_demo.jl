# === 5-ALGORITHM SNF DEMO WITH TIMING ===
# =========================================

println("="^70)
println("SMITH NORMAL FORM - 5 ALGORITHM COMPARISON")
println("With Execution Time")
println("="^70)

# Include all algorithms
include("src/snf_eea.jl")
include("src/snf_divide_conquer.jl")
include("src/snf_modular.jl")
include("src/snf_fraction_free.jl")
include("src/snf_determinant.jl")

using .SmithNormalForm
using .DivideConquerSNF
using .ModularSNF
using .FractionFreeSNF
using .DeterminantBasedSNF

# WARMUP - compile all functions first
println("\n[Warming up Julia JIT compiler...]")
warmup = [1 2; 3 4]
smith_normal_form(warmup)
snf_divide_conquer(warmup)
modular_snf(warmup)
fraction_free_snf(warmup)
determinant_based_snf(warmup)
println("[Done - now measuring actual computation times]\n")

# Test matrices - fixed values for reproducibility
test_matrices = [
    ("2x2 Classic", [6 8; 10 14]),
    ("3×3 Classic", [2 4 4; -6 6 12; 10 4 16]),
    ("4×4 Square", [6 8 14 10; 4 10 4 2; 2 6 12 8; 8 12 6 4]),
    ("5×5 Prime", [2 3 5 7 11; 3 5 7 11 13; 5 7 11 13 17; 7 11 13 17 19; 11 13 17 19 23]),
    ("6×6 Fixed", [3 -7 2 5 -1 8; -4 6 -9 1 7 -2; 5 -3 8 -6 2 4; -1 9 -4 7 -5 3; 6 -2 5 -8 4 1; -7 4 -1 9 -3 6]),
    ("7×7 Fixed", [2 -5 3 7 -1 4 -6; -3 8 -2 5 6 -7 1; 4 -1 9 -3 2 5 -8; -6 7 -4 1 8 -2 3; 5 -9 2 -6 3 7 -4; -1 4 -8 2 -5 9 6; 7 -3 6 -1 4 -8 2]),
    ("8×8 Fixed", [1 -3 5 -7 2 4 -6 8; -2 4 -6 8 -1 3 5 -7; 3 -5 7 -1 4 -6 8 2; -4 6 -8 2 -3 5 7 -1; 5 -7 1 -3 6 -8 2 4; -6 8 -2 4 -5 7 1 -3; 7 -1 3 -5 8 -2 4 6; -8 2 -4 6 -7 1 3 -5])
]

algorithms = [
    ("Classical EEA", (A) -> smith_normal_form(A)),
    ("Divide-Conquer", (A) -> snf_divide_conquer(A)),
    ("Modular", (A) -> modular_snf(A)),
    ("Fraction-Free", (A) -> fraction_free_snf(A)),
    ("Determinant", (A) -> determinant_based_snf(A))
]

for (mat_name, A) in test_matrices
    m, n = size(A)
    min_dim = min(m, n)
    
    println("="^70)
    println("TEST: $mat_name")
    println("="^70)
    
    S_ref, _, _ = smith_normal_form(A)
    eds_ref = [S_ref[i,i] for i in 1:min_dim]
    println("Size: $(m)×$(n) | Expected: ", eds_ref)
    
    println("┌────────────────┬────────────────────────────┬──────────┬────────┐")
    println("│ Algorithm      │ Elementary Divisors        │   Time   │ Status │")
    println("├────────────────┼────────────────────────────┼──────────┼────────┤")
    
    for (name, fn) in algorithms
        try
            t = @elapsed begin
                S, _, _ = fn(A)
            end
            eds = [S[i,i] for i in 1:min(size(S)...)]
            status = eds == eds_ref ? "✅" : "❌"
            time_str = t < 1 ? "$(round(t*1000, digits=1))ms" : "$(round(t, digits=2))s"
            println("│ $(rpad(name,14)) │ $(rpad(string(eds),26)) │ $(rpad(time_str,8)) │   $status   │")
        catch e
            println("│ $(rpad(name,14)) │ Limit: ≤5×5              │ -        │   ⏭️   │")
        end
    end
    println("└────────────────┴────────────────────────────┴──────────┴────────┘\n")
end

println("="^70)
println("ALL TESTS COMPLETED!")
println("="^70)
