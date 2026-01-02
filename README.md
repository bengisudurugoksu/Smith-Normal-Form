# Smith Normal Form Algorithms in Julia

This repository presents a comprehensive study and implementation of **Smith Normal Form (SNF)** algorithms for integer matrices. The project focuses on both classical approaches and modern algorithmic improvements, with an emphasis on exact arithmetic, performance, and scalability.

All algorithms are implemented in **Julia** and evaluated through systematic benchmarking and external validation using **SageMath**.

---

## 📌 Overview

The Smith Normal Form is a canonical diagonal form of integer matrices obtained using unimodular row and column operations. It plays a central role in symbolic computation, algebraic number theory, and module theory, with applications such as:
- Solving linear Diophantine equations  
- Classifying finitely generated abelian groups  
- Analyzing modules over principal ideal domains  

Despite its theoretical simplicity, computing SNF efficiently is challenging due to issues such as **coefficient explosion** in classical algorithms.

---

## 🧠 Implemented Algorithms

This repository includes implementations of the following SNF approaches:

### 1. Classical EEA-Based SNF
- Uses the Extended Euclidean Algorithm (EEA)
- Fast for small matrices
- Suffers from coefficient explosion for larger inputs

### 2. Modular SNF
- Performs computations modulo small primes
- Reconstructs invariant factors using the Chinese Remainder Theorem (CRT)
- Avoids coefficient explosion and scales well to larger matrices

### 3. Fraction-Free SNF
- Eliminates intermediate divisions
- Uses integer arithmetic and gcd-based normalization
- Improves numerical stability and practical performance

### 4. Divide-and-Conquer SNF
- Recursively reduces the matrix size
- Uses dynamic pivot selection
- Provides algorithmic clarity but introduces recursion overhead

### 5. Determinant-Based SNF
- Computes invariant factors directly from minors
- Conceptually simple and definition-based
- Computationally expensive and mainly of theoretical interest

---

## ⚙️ Implementation Details

- Language: **Julia**
- Arithmetic: Exact integer arithmetic
- Design: Modular and readable code structure
- Focus: Correctness, clarity, and algorithmic comparison

Each algorithm follows the mathematical definition of SNF and preserves unimodularity throughout the transformations.

---

## 📊 Benchmarking & Evaluation

The repository includes benchmarking experiments that measure:
- Execution time vs. matrix size
- Practical scalability of each algorithm
- Performance trade-offs between methods

Results show that **modular SNF** provides the best overall performance for medium and large matrices, while classical and fraction-free methods remain useful for small-scale problems.

---

## ✅ Validation with SageMath

To ensure correctness, all implementations were cross-validated using **SageMath**, a widely accepted computer algebra system.  
Identical input matrices were tested in both Julia and SageMath, yielding a **100% match** in invariant factors across all test cases.

---

## 🎯 Key Takeaways

- No single SNF algorithm is optimal for all cases
- Modular SNF is the most scalable and practical choice
- Fraction-free methods offer a strong balance between stability and performance
- Determinant-based methods are best suited for theoretical exploration

---

## 📚 References

Relevant academic references and resources are listed in the accompanying presentation and documentation.

---

## 👩‍💻 Authors

This project was developed as part of an academic study on Smith Normal Form algorithms and their computational properties.
