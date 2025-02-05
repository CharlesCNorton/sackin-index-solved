"""
Rigorous Testing of the Sackin Index Closed–Form Generating Function

This script performs the following:
1. Computes the total Sackin index S(n) for full binary trees with n leaves 
   using a dynamic programming approach based on the recurrence:
       S(1) = 0,
       S(n) = sum_{i=1}^{n-1}[ S(i)*T(n-i) + S(n-i)*T(i) + n*T(i)*T(n-i) ]
   where T(n) is the number of full binary trees with n leaves 
   (with T(1)=1 and, for n>=2, T(n) equals the (n-1)th Catalan number).

2. Computes T(n) using the well-known formula for Catalan numbers:
       T(n) = 1 for n=1, and for n>=2, T(n) = Catalan(n-1) = binomial(2n-2, n-1)/(n)

3. Extracts coefficients from the closed–form generating function:
       F_S(x) = x*(1 - sqrt(1-4*x))/(1-4*x)
   using a series expansion (up to n = 50).

4. Compares the two sets of S(n) values and prints a table.

5. Plots the average Sackin index per tree (i.e. S(n)/T(n)) on a log–log plot to confirm
   the expected asymptotic behavior.
"""

import sympy as sp
import math
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Setup: Define symbols and initialize printing
# -----------------------------
sp.init_printing(use_unicode=True)
x = sp.symbols('x', positive=True, real=True)

# -----------------------------
# Part 1: Define Catalan Number and T(n)
# -----------------------------
def catalan(n):
    """
    Compute the nth Catalan number.
    Note: For full binary trees with n leaves, we need Catalan(n-1) for n>=2,
    with the convention T(1)=1.
    """
    return sp.binomial(2*n, n) // (n+1)

def T_n(n):
    """
    Return the number of full binary trees with n leaves.
    T(1)=1; for n>=2, T(n) = Catalan(n-1)
    """
    if n == 1:
        return sp.Integer(1)
    else:
        return catalan(n-1)

# -----------------------------
# Part 2: Compute S(n) via Dynamic Programming (Recurrence)
# -----------------------------
def compute_Sackin_index(N_max):
    """
    Computes the total Sackin index S(n) for n=1 to N_max using the recurrence:
      S(1) = 0,
      S(n) = sum_{i=1}^{n-1}[ S(i)*T(n-i) + S(n-i)*T(i) + n*T(i)*T(n-i) ]
    Returns:
      S_dict: dictionary mapping n -> S(n)
      T_dict: dictionary mapping n -> T(n)
    """
    S_dict = {1: sp.Integer(0)}
    T_dict = {1: sp.Integer(1)}
    # Precompute T(n) for n=1 to N_max
    for n in range(2, N_max+1):
        T_dict[n] = T_n(n)
    # Compute S(n) for n>=2
    for n in range(2, N_max+1):
        s_val = sp.Integer(0)
        for i in range(1, n):
            # Contribution from left subtree (of size i) and right subtree (of size n-i)
            term1 = S_dict[i] * T_dict[n-i]
            term2 = S_dict[n-i] * T_dict[i]
            # Every leaf increases its depth by 1 when subtrees are combined: contributes n*T(i)*T(n-i)
            term3 = n * T_dict[i] * T_dict[n-i]
            s_val += term1 + term2 + term3
        S_dict[n] = sp.simplify(s_val)
    return S_dict, T_dict

# -----------------------------
# Part 3: Extract Coefficients from the Closed-Form Generating Function
# -----------------------------
def series_coefficients_F_S(N_max):
    """
    Expand the closed-form generating function:
         F_S(x) = x*(1 - sqrt(1-4*x))/(1-4*x)
    and extract coefficients for x^n, for n = 1 to N_max.
    Returns a dictionary mapping n -> coefficient.
    """
    F_S = x*(1 - sp.sqrt(1-4*x))/(1-4*x)
    # Expand the series to the required order.
    # We request series up to order N_max+1 to get coefficients for x^1,...,x^N_max.
    series_expansion = sp.series(F_S, x, 0, N_max+1).removeO().expand()
    coeff_dict = {}
    for n in range(1, N_max+1):
        coeff = sp.nsimplify(series_expansion.coeff(x, n))
        coeff_dict[n] = coeff
    return coeff_dict

# -----------------------------
# Part 4: Compare the Two Methods and Output the Results
# -----------------------------
def compare_results(N_max):
    """
    Computes S(n) via dynamic programming and via series extraction,
    then compares them for n=1 to N_max.
    Returns a pandas DataFrame with columns:
      n, T(n), S(n) from recurrence, S(n) from generating function, and Delta.
    """
    S_recur, T_dict = compute_Sackin_index(N_max)
    S_series = series_coefficients_F_S(N_max)
    
    data = []
    for n in range(1, N_max+1):
        s_rec = S_recur[n]
        s_ser = S_series[n]
        delta = sp.simplify(s_ser - s_rec)
        data.append({
            "n": n,
            "#Trees T(n)": T_dict[n],
            "S(n) (Recurrence)": s_rec,
            "S(n) (Series Coefficient)": s_ser,
            "Delta": delta
        })
    
    df = pd.DataFrame(data)
    return df

# -----------------------------
# Part 5: Plotting the Average Sackin Index per Tree
# -----------------------------
def plot_average_Sackin(N_max, S_recur, T_dict):
    """
    Plots the average Sackin index per tree, i.e., A(n) = S(n)/T(n), on a log-log plot.
    """
    avg_S = []
    n_vals = list(range(1, N_max+1))
    for n in n_vals:
        avg = sp.N(S_recur[n] / T_dict[n])
        avg_S.append(float(avg))
    
    plt.figure(figsize=(10, 6))
    # Set the scales explicitly (base 10)
    plt.xscale('log', base=10)
    plt.yscale('log', base=10)
    plt.plot(n_vals, avg_S, 'bo-')
    plt.xlabel("n (Number of Leaves)")
    plt.ylabel("Average Sackin Index per Tree, A(n) = S(n)/T(n)")
    plt.title("Log-Log Plot of Average Sackin Index per Tree")
    plt.grid(True, which="both", ls="--")
    plt.show()

# -----------------------------
# Main Testing Routine: Test up to n = 50
# -----------------------------
def main():
    N_max = 50  # Testing up to 50 leaves
    # Compute values via dynamic programming and series extraction
    df_results = compare_results(N_max)
    print("Comparison of S(n) from Recurrence and from Series Expansion:")
    print(df_results.to_string(index=False))
    
    # Check if all Delta values are 0 (i.e., perfect match)
    discrepancies = df_results[df_results["Delta"] != 0]
    if discrepancies.empty:
        print("\nAll coefficients match exactly between the two methods up to n =", N_max)
    else:
        print("\nDiscrepancies found:")
        print(discrepancies)
    
    # Retrieve S(n) and T(n) for plotting the average Sackin index per tree.
    S_recur, T_dict = compute_Sackin_index(N_max)
    plot_average_Sackin(N_max, S_recur, T_dict)
    
if __name__ == "__main__":
    main()
