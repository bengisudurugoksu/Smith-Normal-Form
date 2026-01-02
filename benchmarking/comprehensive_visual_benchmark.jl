# === COMPREHENSIVE VISUAL BENCHMARK ===
# Creates tables and visualizations for presentation
# ==========================================

using Plots, Printf, Statistics
using LinearAlgebra

println("="^80)
println("COMPREHENSIVE SNF BENCHMARK WITH VISUALIZATIONS")
println("="^80)

# Include algorithms
include("src/snf_eea.jl")
include("src/snf_divide_conquer.jl")
include("src/snf_modular_nemo.jl")
include("src/snf_fraction_free.jl")
include("src/snf_determinant.jl")

using .SmithNormalForm
using .DivideConquerSNF
using .ModularNemoSNF
using .FractionFreeSNF
using .DeterminantBasedSNF

# Warmup
println("\n[Warming up...]")
warmup = [1 2; 3 4]
smith_normal_form(warmup)
snf_divide_conquer(warmup)
snf_modular_nemo(warmup)
fraction_free_snf(warmup)
determinant_based_snf(warmup)
println("[Done]\n")

# Test matrices with varied sizes
test_matrices = [
    ("3×3", [2 4 4; -6 6 12; 10 4 16]),
    ("4×4", [6 8 14 10; 4 10 4 2; 2 6 12 8; 8 12 6 4]),
    ("5×5", [2 3 5 7 11; 3 5 7 11 13; 5 7 11 13 17; 7 11 13 17 19; 11 13 17 19 23]),
    ("6×6", [3 -7 2 5 -1 8; -4 6 -9 1 7 -2; 5 -3 8 -6 2 4; -1 9 -4 7 -5 3; 6 -2 5 -8 4 1; -7 4 -1 9 -3 6]),
    ("7×7", [2 -5 3 7 -1 4 -6; -3 8 -2 5 6 -7 1; 4 -1 9 -3 2 5 -8; -6 7 -4 1 8 -2 3; 5 -9 2 -6 3 7 -4; -1 4 -8 2 -5 9 6; 7 -3 6 -1 4 -8 2]),
    ("8×8", [1 -3 5 -7 2 4 -6 8; -2 4 -6 8 -1 3 5 -7; 3 -5 7 -1 4 -6 8 2; -4 6 -8 2 -3 5 7 -1; 5 -7 1 -3 6 -8 2 4; -6 8 -2 4 -5 7 1 -3; 7 -1 3 -5 8 -2 4 6; -8 2 -4 6 -7 1 3 -5]),
    ("10×10", rand(-10:10, 10, 10))
]

algorithms = [
    ("Classical EEA", smith_normal_form, "blue"),
    ("Divide-Conquer", snf_divide_conquer, "red"),
    ("Modular (Nemo)", snf_modular_nemo, "green"),
    ("Fraction-Free", fraction_free_snf, "purple"),
    ("Determinant", determinant_based_snf, "orange")
]


# Storage for results
results = Dict()
sage_results = Dict()
memory_results = Dict()

println("="^80)
println("PHASE 1: Running Julia Algorithms")
println("="^80)

for (mat_name, A) in test_matrices
    println("\n📊 Testing $mat_name matrix...")
    results[mat_name] = Dict()
    memory_results[mat_name] = Dict()
    
    for (name, fn, color) in algorithms
        try
            # Time measurement
            times = Float64[]
            for _ in 1:5  # Run 5 times for stability
                t = @elapsed fn(A)
                push!(times, t)
            end
            avg_time = mean(times)
            std_time = std(times)
            
            # Memory measurement
            mem_stats = @timed fn(A)
            mem_allocs = mem_stats.bytes
            
            results[mat_name][name] = (avg_time, std_time)
            memory_results[mat_name][name] = mem_allocs
            
            println("   ✓ $name: $(round(avg_time*1000, digits=2))ms (±$(round(std_time*1000, digits=2))ms)")
        catch e
            println("   ✗ $name: FAILED")
            results[mat_name][name] = (NaN, NaN)
            memory_results[mat_name][name] = 0
        end
    end
end

println("\n" * "="^80)
println("PHASE 2: Running SageMath Comparison")
println("="^80)

for (mat_name, A) in test_matrices
    println("\n📊 Testing $mat_name matrix with SageMath...")
    m, n = size(A)
    
    # Convert matrix to SageMath format
    mat_str = "["
    for i in 1:m
        mat_str *= "["
        for j in 1:n
            mat_str *= string(A[i,j])
            j < n && (mat_str *= ", ")
        end
        mat_str *= "]"
        i < m && (mat_str *= ", ")
    end
    mat_str *= "]"
    
    sage_cmd = """
    import time
    times = []
    for _ in range(5):
        start = time.time()
        M = matrix(ZZ, $mat_str)
        eds = M.elementary_divisors()
        elapsed = time.time() - start
        times.append(elapsed)
    
    import statistics
    avg_time = statistics.mean(times)
    std_time = statistics.stdev(times) if len(times) > 1 else 0
    print(f"TIME:{avg_time}")
    print(f"STD:{std_time}")
    """
    
    # Save and run
    open("/tmp/sage_bench.sage", "w") do f
        write(f, sage_cmd)
    end
    
    try
        sage_output = read(`sage /tmp/sage_bench.sage`, String)
        
        sage_time = NaN
        sage_std = NaN
        for line in split(sage_output, "\n")
            if startswith(line, "TIME:")
                sage_time = parse(Float64, replace(line, "TIME:" => ""))
            elseif startswith(line, "STD:")
                sage_std = parse(Float64, replace(line, "STD:" => ""))
            end
        end
        
        sage_results[mat_name] = (sage_time, sage_std)
        println("   ✓ SageMath: $(round(sage_time*1000, digits=2))ms (±$(round(sage_std*1000, digits=2))ms)")
    catch e
        println("   ✗ SageMath: FAILED")
        sage_results[mat_name] = (NaN, NaN)
    end
end

println("\n" * "="^80)
println("PHASE 3: Generating Visualizations")
println("="^80)

# Extract data for plotting
matrix_sizes = [parse(Int, split(name, "×")[1]) for (name, _) in test_matrices]
alg_names = [name for (name, _, _) in algorithms]

# 1. Time Comparison Plot
println("\n📈 Creating time comparison plot...")
p1 = plot(title="Execution Time vs Matrix Size", xlabel="Matrix Size (n×n)", 
         ylabel="Time (ms)", yscale=:log10, legend=:topleft, size=(900, 600),
         dpi=300)

for (name, fn, color) in algorithms
    times = [get(get(results, mat_name, Dict()), name, (NaN, NaN))[1] * 1000 
             for (mat_name, _) in test_matrices]
    plot!(p1, matrix_sizes, times, label=name, marker=:circle, 
          linewidth=2, markersize=6, color=Symbol(color))
end

# Add SageMath
sage_times = [get(sage_results, mat_name, (NaN, NaN))[1] * 1000 
              for (mat_name, _) in test_matrices]
plot!(p1, matrix_sizes, sage_times, label="SageMath", marker=:diamond, 
      linewidth=2, markersize=6, color=:black, linestyle=:dash)

savefig(p1, "benchmark_time_comparison.png")
println("   ✓ Saved: benchmark_time_comparison.png")

# 2. Speedup vs SageMath
println("\n📈 Creating speedup comparison plot...")
p2 = plot(title="Julia Speedup vs SageMath", xlabel="Matrix Size (n×n)", 
         ylabel="Speedup Factor (higher is better)", legend=:topright, 
         size=(900, 600), dpi=300)

for (name, fn, color) in algorithms
    speedups = Float64[]
    for (mat_name, _) in test_matrices
        julia_time = get(get(results, mat_name, Dict()), name, (NaN, NaN))[1]
        sage_time = get(sage_results, mat_name, (NaN, NaN))[1]
        if !isnan(julia_time) && !isnan(sage_time) && julia_time > 0
            push!(speedups, sage_time / julia_time)
        else
            push!(speedups, NaN)
        end
    end
    plot!(p2, matrix_sizes, speedups, label=name, marker=:circle, 
          linewidth=2, markersize=6, color=Symbol(color))
end

hline!(p2, [1.0], label="Equal Performance", color=:black, linestyle=:dot, linewidth=2)
savefig(p2, "benchmark_speedup_comparison.png")
println("   ✓ Saved: benchmark_speedup_comparison.png")

# 3. Memory Usage Comparison
println("\n📈 Creating memory usage plot...")
p3 = plot(title="Memory Allocations vs Matrix Size", xlabel="Matrix Size (n×n)", 
         ylabel="Memory (MB)", legend=:topleft, size=(900, 600), dpi=300)

for (name, fn, color) in algorithms
    mem_vals = [get(get(memory_results, mat_name, Dict()), name, 0) / 1024^2 
                for (mat_name, _) in test_matrices]
    plot!(p3, matrix_sizes, mem_vals, label=name, marker=:circle, 
          linewidth=2, markersize=6, color=Symbol(color))
end

savefig(p3, "benchmark_memory_comparison.png")
println("   ✓ Saved: benchmark_memory_comparison.png")

# 4. Bar Chart - Average Performance Across All Tests
println("\n📈 Creating average performance bar chart...")
avg_times = Float64[]
for (name, fn, color) in algorithms
    times = [get(get(results, mat_name, Dict()), name, (NaN, NaN))[1] * 1000 
             for (mat_name, _) in test_matrices]
    push!(avg_times, mean(filter(!isnan, times)))
end

sage_avg = mean(filter(!isnan, [get(sage_results, mat_name, (NaN, NaN))[1] * 1000 
                                  for (mat_name, _) in test_matrices]))

p4 = bar([alg_names; "SageMath"], [avg_times; sage_avg], 
         title="Average Execution Time", ylabel="Time (ms)", 
         xlabel="Algorithm", legend=false, color=[:blue :red :green :purple :orange :black],
         size=(900, 600), dpi=300)

savefig(p4, "benchmark_average_comparison.png")
println("   ✓ Saved: benchmark_average_comparison.png")

# 5. Generate LaTeX Table
println("\n📄 Creating LaTeX table...")
global latex_output = """
\\begin{table}[h]
\\centering
\\caption{Julia vs SageMath - Execution Time Comparison (ms)}
\\begin{tabular}{|l|c|c|c|c|c|c|}
\\hline
\\textbf{Size} & \\textbf{Classical} & \\textbf{Divide-Conq} & \\textbf{Modular} & \\textbf{Fraction-Free} & \\textbf{Determinant} & \\textbf{SageMath} \\\\
\\hline
"""

for (mat_name, _) in test_matrices
    global latex_output
    row = mat_name * " & "
    for (name, _, _) in algorithms
        time_val = get(get(results, mat_name, Dict()), name, (NaN, NaN))[1] * 1000
        row *= @sprintf("%.2f", time_val) * " & "
    end
    sage_time = get(sage_results, mat_name, (NaN, NaN))[1] * 1000
    row *= @sprintf("%.2f", sage_time) * " \\\\\n"
    latex_output *= row
end

latex_output *= """
\\hline
\\end{tabular}
\\end{table}
"""

open("benchmark_table.tex", "w") do f
    write(f, latex_output)
end
println("   ✓ Saved: benchmark_table.tex")

# 6. Generate CSV for Excel/Sheets
println("\n📊 Creating CSV file...")
open("benchmark_results.csv", "w") do f
    # Header
    write(f, "Matrix Size,")
    for (name, _, _) in algorithms
        write(f, "$name Time (ms),$name Memory (MB),")
    end
    write(f, "SageMath Time (ms),Best Julia,Speedup\n")
    
    # Data rows
    for (mat_name, _) in test_matrices
        write(f, "$mat_name,")
        
        best_time = Inf
        best_alg = ""
        
        for (name, _, _) in algorithms
            time_val = get(get(results, mat_name, Dict()), name, (NaN, NaN))[1] * 1000
            mem_val = get(get(memory_results, mat_name, Dict()), name, 0) / 1024^2
            write(f, @sprintf("%.3f,%.2f,", time_val, mem_val))
            
            if time_val < best_time
                best_time = time_val
                best_alg = name
            end
        end
        
        sage_time = get(sage_results, mat_name, (NaN, NaN))[1] * 1000
        speedup = sage_time / best_time
        write(f, @sprintf("%.3f,%s,%.2fx\n", sage_time, best_alg, speedup))
    end
end
println("   ✓ Saved: benchmark_results.csv")

# Print Summary Table for Console
println("\n" * "="^80)
println("RESULTS SUMMARY")
println("="^80)

println("\n┌─────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐")
println("│  Size   │  Classical   │ Divide-Conq  │   Modular    │ Fraction-Free│ Determinant  │   SageMath   │")
println("├─────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤")

for (mat_name, _) in test_matrices
    print("│ $(rpad(mat_name, 7)) │")
    
    for (name, _, _) in algorithms
        time_val = get(get(results, mat_name, Dict()), name, (NaN, NaN))[1] * 1000
        time_str = @sprintf("%.2f ms", time_val)
        print(" $(rpad(time_str, 12)) │")
    end
    
    sage_time = get(sage_results, mat_name, (NaN, NaN))[1] * 1000
    sage_str = @sprintf("%.2f ms", sage_time)
    println(" $(rpad(sage_str, 12)) │")
end

println("└─────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘")

# Performance Summary
println("\n" * "="^80)
println("PERFORMANCE ANALYSIS")
println("="^80)

total_wins = Dict(name => 0 for (name, _, _) in algorithms)

for (mat_name, _) in test_matrices
    best_time = Inf
    best_alg = ""
    
    for (name, _, _) in algorithms
        time_val = get(get(results, mat_name, Dict()), name, (NaN, NaN))[1]
        if !isnan(time_val) && time_val < best_time
            best_time = time_val
            best_alg = name
        end
    end
    
    if best_alg != ""
        total_wins[best_alg] += 1
        sage_time = get(sage_results, mat_name, (NaN, NaN))[1]
        speedup = sage_time / best_time
        println("📊 $mat_name: $best_alg won ($(round(speedup, digits=2))× faster than SageMath)")
    end
end

println("\n🏆 ALGORITHM WINS:")
for (name, wins) in sort(collect(total_wins), by=x->x[2], rev=true)
    println("   $name: $wins tests")
end

println("\n" * "="^80)
println("✅ BENCHMARK COMPLETED!")
println("="^80)
println("\nGenerated files:")
println("  • benchmark_time_comparison.png")
println("  • benchmark_speedup_comparison.png")
println("  • benchmark_memory_comparison.png")
println("  • benchmark_average_comparison.png")
println("  • benchmark_table.tex (for LaTeX)")
println("  • benchmark_results.csv (for Excel/Sheets)")
