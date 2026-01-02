using LinearAlgebra

"""
    smith_normal_form_with_uv(M)

Computes the Smith Normal Form (SNF) of an integer matrix M using a 
Divide and Conquer (Recursive) approach.

Returns: (S, U, V) where S = U * M * V
- S: The diagonal Smith Normal Form matrix.
- U: The left unimodular transformation matrix.
- V: The right unimodular transformation matrix.

Reference: Project Proposal, Section "Algorithm for Computing the SNF"[cite: 18].
"""
function smith_normal_form_with_uv(M::Matrix{T}) where T <: Integer
    # Convert input to BigInt to prevent integer overflow during intermediate calculations.
    # As noted in the proposal, coefficient growth can occur[cite: 25].
    A = BigInt.(copy(M))
    m, n = size(A)
    
    # Initialize transformation matrices U and V as Identity matrices.
    U = Matrix{BigInt}(I, m, m)
    V = Matrix{BigInt}(I, n, n)
    
    # --- BASE CASE ---
    # If the matrix dimension is 0, the recursion terminates.
    if m == 0 || n == 0
        return A, U, V
    end

    # --- HELPER FUNCTION: PIVOT SELECTION ---
    # Finds the entry with the smallest non-zero absolute value in the matrix
    # and moves it to the (1,1) position using swaps[cite: 20].
    function move_best_pivot_to_top!(Mat, MatU, MatV, rows, cols)
        min_val = BigInt(0)
        p_r, p_c = 0, 0
        found = false
        
        # Search the entire matrix for the best pivot
        for i in 1:rows, j in 1:cols
            val = abs(Mat[i, j])
            if val != 0
                if !found || val < min_val
                    min_val = val
                    p_r, p_c = i, j
                    found = true
                end
            end
        end
        
        if !found return false end # Matrix is entirely zero
        
        # Swap Rows (affects A and U)
        if p_r != 1
            Mat[1, :], Mat[p_r, :] = Mat[p_r, :], Mat[1, :]
            MatU[1, :], MatU[p_r, :] = MatU[p_r, :], MatU[1, :]
        end
        
        # Swap Columns (affects A and V)
        if p_c != 1
            Mat[:, 1], Mat[:, p_c] = Mat[:, p_c], Mat[:, 1]
            MatV[:, 1], MatV[:, p_c] = MatV[:, p_c], MatV[:, 1]
        end
        return true
    end

    # --- STEP 1: ELIMINATION LOOP (REDUCTION) ---
    # Iteratively clear the first row and first column using Euclidean steps[cite: 21, 22].
    while true
        # Select best pivot and move to (1,1)
        if !move_best_pivot_to_top!(A, U, V, m, n)
            return A, U, V # Matrix is zero, done.
        end
        
        changed = false
        
        # Column Clearing: Eliminate entries A[2:end, 1]
        # Operation: Row_i = Row_i - q * Row_1
        for i in 2:m
            if A[i, 1] != 0
                q = div(A[i, 1], A[1, 1])
                A[i, :] .-= q .* A[1, :]
                U[i, :] .-= q .* U[1, :] # Apply same operation to U
                if A[i, 1] != 0; changed = true; end
            end
        end
        
        # Row Clearing: Eliminate entries A[1, 2:end]
        # Operation: Col_j = Col_j - q * Col_1
        for j in 2:n
            if A[1, j] != 0
                q = div(A[1, j], A[1, 1])
                A[:, j] .-= q .* A[:, 1]
                V[:, j] .-= q .* V[:, 1] # Apply same operation to V
                if A[1, j] != 0; changed = true; end
            end
        end
        
        # If no changes occurred, the pivot A[1,1] divides all entries in its row/col.
        if !changed
            break
        end
    end

    # --- STEP 2: RECURSION (DIVIDE & CONQUER) ---
    # Isolate the submatrix (bottom-right block) and solve recursively[cite: 23].
    sub_A = A[2:end, 2:end]
    sub_S, sub_U, sub_V = smith_normal_form_with_uv(sub_A)
    
    # --- STEP 3: MERGE RESULTS ---
    # Embed the recursive results back into the main matrices.
    
    # Update A with the diagonalized submatrix
    A[2:end, 2:end] = sub_S
    
    # Update global U: Combine current operations with sub-problem operations
    U_lift = Matrix{BigInt}(I, m, m)
    U_lift[2:end, 2:end] = sub_U
    U = U_lift * U
    
    # Update global V: Combine current operations with sub-problem operations
    V_lift = Matrix{BigInt}(I, n, n)
    V_lift[2:end, 2:end] = sub_V
    V = V * V_lift

    # --- STEP 4: DIVISIBILITY CHECK (INVARIANT FACTORS) ---
    # Ensure d_1 divides d_2 (A[1,1] | A[2,2])[cite: 23].
    if min(m, n) > 1
        d1 = A[1, 1]
        d2 = A[2, 2]
        if d2 % d1 != 0
            # If divisibility fails, add row 2 to row 1 to create a GCD interaction
            A[1, :] .+= A[2, :]
            U[1, :] .+= U[2, :] 
            
            # Restart the process for the modified matrix (Tail Recursion logic)
            res_S, res_U, res_V = smith_normal_form_with_uv(A)
            
            # Merge the new transformation matrices
            return res_S, res_U * U, V * res_V
        end
    end

    # --- FINALIZATION: SIGN CORRECTION ---
    # Ensure the invariant factor is positive (unique up to units)[cite: 11].
    if A[1, 1] < 0
        A[1, :] .*= -1
        U[1, :] .*= -1
    end

    return A, U, V
end