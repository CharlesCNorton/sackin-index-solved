"""
This script rigorously verifies the closed-form generating function for the Sackin index of full binary trees.
It computes the Sackin index S(n) for n = 1, 2, …, N in two ways:

  1. Using a recurrence (dynamic programming):
       S(1) = 0,
       For n ≥ 2:
         S(n) = Σ_{i=1}^{n-1} [ S(i) * T(n-i) + S(n-i) * T(i) + n * T(i) * T(n-i) ]
     where T(n) is the number of full binary trees with n leaves:
       T(1) = 1, and for n ≥ 2, T(n) = Catalan(n-1).

  2. By expanding the candidate generating function:
         F_S(x) = x*(1 - sqrt(1-4*x))/(1-4*x)
     whose series expansion should yield coefficients S(n).

The script then compares the two computations for each n.
For large N the computation may be computationally intensive.

Dependencies:
  - sympy
  - pandas
  - (Optional) time, sys

Run this script from the command line and enter a positive integer when prompted.
"""

import sympy as sp
import pandas as pd
import sys
import time

# Initialize sympy printing and define symbol x.
sp.init_printing(use_unicode=True)
x = sp.symbols('x', real=True, positive=True)

# Define the generating function for full binary trees (shifted Catalan numbers):
# T(x) = (1 - sqrt(1-4*x))/2.
T_gen = (1 - sp.sqrt(1-4*x)) / 2

# Candidate generating function for the Sackin index:
F_S_candidate = x * (1 - sp.sqrt(1-4*x)) / (1-4*x)

def catalan(n):
    """
    Compute the nth Catalan number:
        Catalan(n) = binomial(2*n, n) / (n+1).
    """
    return sp.binomial(2 * n, n) // (n + 1)

def T_n(n):
    """
    Return the number of full binary trees with n leaves.
    T(1) = 1, and for n >= 2, T(n) = Catalan(n-1).
    """
    if n == 1:
        return sp.Integer(1)
    else:
        return catalan(n - 1)

def compute_S_n(max_n):
    """
    Compute S(n) for n = 1,2,...,max_n using the recurrence:
      S(1) = 0,
      S(n) = Σ_{i=1}^{n-1} [ S(i)*T(n-i) + S(n-i)*T(i) + n * T(i)*T(n-i) ].
    Returns a dictionary mapping n -> S(n).
    """
    S = {1: sp.Integer(0)}
    for n in range(2, max_n + 1):
        s_val = sp.Integer(0)
        for i in range(1, n):
            s_val += S[i] * T_n(n - i) + S[n - i] * T_n(i) + n * T_n(i) * T_n(n - i)
        S[n] = sp.simplify(s_val)
        print(f"Computed S({n})")
    return S

def compute_candidate_series_coeffs(max_n):
    """
    Compute the coefficients for F_S(x) (the candidate generating function for the Sackin index)
    from x^1 up to x^(max_n). Returns a list of coefficients.
    """
    order = max_n + 1  # we need terms up to x^(max_n)
    series_expansion = sp.series(F_S_candidate, x, 0, order).removeO().expand()
    coeffs = []
    for n in range(1, max_n + 1):
        coeff = sp.nsimplify(series_expansion.coeff(x, n))
        coeffs.append(coeff)
        print(f"Candidate coefficient for x^{n} computed.")
    return coeffs

def main():
    try:
        max_n_input = input("Enter the maximum value of n to verify (positive integer): ")
        max_n = int(max_n_input)
        if max_n < 1:
            print("Please enter an integer greater than or equal to 1.")
            sys.exit(1)
    except Exception as e:
        print("Invalid input. Please enter a positive integer.")
        sys.exit(1)

    start_time = time.time()
    print("\nStarting computation of S(n) via dynamic programming...")
    S_values = compute_S_n(max_n)
    print("Dynamic programming computation completed.\n")

    print("Computing series coefficients from the candidate generating function...")
    candidate_coeffs = compute_candidate_series_coeffs(max_n)
    print("Candidate series coefficients computation completed.\n")

    # Prepare the verification table.
    data = []
    for n in range(1, max_n + 1):
        trees_count = T_n(n)
        S_recur = S_values[n]
        S_cand = candidate_coeffs[n - 1]
        delta = sp.simplify(S_cand - S_recur)
        data.append({
            'n': n,
            '#Trees': trees_count,
            'Sackin (Recurrence)': S_recur,
            'Sackin (Candidate)': S_cand,
            'Delta': delta
        })

    df = pd.DataFrame(data)
    print("Verification Results:")
    print(df.to_string(index=False))

    # Check if all differences (Delta) are zero.
    if all(delta == 0 for delta in df['Delta']):
        print(f"\nSuccess: The candidate generating function exactly matches the recurrence-based values for all n up to {max_n}.")
    else:
        print("\nError: Discrepancies found between the candidate generating function and the recurrence-based values.")

    elapsed = time.time() - start_time
    print(f"\nTotal computation time: {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()
