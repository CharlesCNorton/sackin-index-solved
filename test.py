"""
Torture Testing the Sackin Index Solution

This script performs a rigorous verification of the closed-form generating function
for the Sackin index of full binary trees. We verify our candidate generating function

    F_S(x) = x*(1-sqrt(1-4x))/(1-4x)

by comparing its series coefficients with values computed from a recurrence.

The recurrence is derived as follows:
  - Let T(n) denote the number of full binary trees with n leaves 
    (T(1)=1 and for n>=2, T(n) is the (n-1)th Catalan number).
  - Let S(n) be the total Sackin index (sum of depths of all leaves, with root at depth 0)
    summed over all trees with n leaves.
  - For n=1, S(1)=0.
  - For n>=2, a full binary tree is formed by joining two subtrees with i and n-i leaves.
    Each such joining increases the depth of every leaf by 1 (contributing n overall).
    Hence, the recurrence is:
    
      S(n) = Sum_{i=1}^{n-1} [ S(i)*T(n-i) + S(n-i)*T(i) + n * T(i)*T(n-i) ].

We compute S(n) for n up to N_max using dynamic programming and compare with
the series coefficients from F_S(x).
"""

import sympy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set up sympy symbols and pretty printing
sp.init_printing(use_unicode=True)
x = sp.symbols('x', real=True, positive=True)

# --------------------------------------------------
# Part 1: Candidate Generating Function for Sackin Index
# --------------------------------------------------
# Our candidate is:
F_S_candidate = x * (1 - sp.sqrt(1 - 4*x)) / (1 - 4*x)

# Expand the candidate series up to a high order.
order_FS = 16  # test for n=1,...,15
F_S_series = sp.series(F_S_candidate, x, 0, order_FS).removeO().expand()
Sackin_candidate_coeffs = [sp.nsimplify(F_S_series.coeff(x, n)) for n in range(1, order_FS)]

print("Candidate Generating Function for Sackin Index (series coefficients):")
for n, coeff in enumerate(Sackin_candidate_coeffs, start=1):
    print("  n = {:2d}: {}".format(n, coeff))

# --------------------------------------------------
# Part 2: Recurrence–Based Computation of Sackin Index Totals
# --------------------------------------------------
def catalan(n):
    """Return the nth Catalan number."""
    return sp.binomial(2*n, n) // (n+1)

def T(n):
    """
    Number of full binary trees with n leaves.
    For n = 1, T(1) = 1, and for n>=2, T(n) = Catalan(n-1).
    """
    if n == 1:
        return sp.Integer(1)
    else:
        return catalan(n-1)

# Use dynamic programming to compute S(n) for n = 1,..., N_max
N_max = 15  # test for n up to 15
S = {1: sp.Integer(0)}  # S(1) = 0

for n in range(2, N_max+1):
    s_val = sp.Integer(0)
    # Sum over splits: for each i from 1 to n-1, the tree is built from a left subtree of size i and right subtree of size n-i.
    for i in range(1, n):
        # The total Sackin index from left subtree: S(i)
        # The total Sackin index from right subtree: S(n-i)
        # Additionally, every leaf in both subtrees gets an extra depth of 1 (n leaves total).
        s_val += S[i] * T(n-i) + S[n-i] * T(i) + n * T(i) * T(n-i)
    S[n] = sp.simplify(s_val)

# Prepare lists for comparison
Sackin_recurrence = [S[n] for n in range(1, N_max+1)]
Trees_counts = [T(n) for n in range(1, N_max+1)]

print("\nRecurrence-Based Sackin Index Totals:")
for n in range(1, N_max+1):
    print("  n = {:2d}: #Trees = {:4d}, S(n) = {}".format(n, int(T(n)), Sackin_recurrence[n-1]))

# --------------------------------------------------
# Part 3: Compare Candidate Coefficients with Recurrence Values
# --------------------------------------------------
n_compare = min(len(Sackin_candidate_coeffs), N_max)
Delta_Sackin = []
for n in range(1, n_compare+1):
    candidate_val = Sackin_candidate_coeffs[n-1]
    recurrence_val = sp.Integer(Sackin_recurrence[n-1])
    delta = sp.simplify(candidate_val - recurrence_val)
    Delta_Sackin.append(delta)

df_sackin = pd.DataFrame({
    "n": list(range(1, n_compare+1)),
    "#Trees": [int(T(n)) for n in range(1, n_compare+1)],
    "Sackin (Recurrence)": [int(Sackin_recurrence[n-1]) for n in range(1, n_compare+1)],
    "Sackin (Candidate)": [candidate for candidate in Sackin_candidate_coeffs[:n_compare]],
    "Delta": [delta for delta in Delta_Sackin]
})
print("\nComparison for Sackin Index:")
print(df_sackin.to_string(index=False))

# Check if any discrepancy exists:
if any(delta != 0 for delta in Delta_Sackin):
    print("\nDiscrepancies detected! Candidate does not exactly match recurrence-based values.")
else:
    print("\nNo discrepancy detected: the candidate generating function matches the recurrence-based values perfectly.")

# --------------------------------------------------
# Part 4: Asymptotic Verification (Optional)
# --------------------------------------------------
# We now plot the average Sackin index (S(n) / T(n)) versus n on a log-log scale.
avg_Sackin = [float(Sackin_recurrence[n-1] / T(n)) for n in range(1, N_max+1)]
n_values = np.array(range(1, N_max+1))

plt.figure(figsize=(8, 5))
plt.loglog(n_values, avg_Sackin, 'bo-', label="Average Sackin (recurrence)")
plt.xlabel("n (number of leaves)")
plt.ylabel("Average Sackin Index")
plt.title("Average Sackin Index vs n (log-log)")
plt.legend()
plt.grid(True)
plt.show()
