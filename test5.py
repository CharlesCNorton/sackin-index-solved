"""
Self-contained script to test the closed-form generating function for the Sackin index
of full binary trees using multiple independent convergence measures.
"""

import sympy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Part 1: Define the Generating Function and Expand Series
# -----------------------------
sp.init_printing(use_unicode=True)
x = sp.symbols('x', real=True, positive=True)

# Closed-form generating function for the Sackin index:
# F_S(x) = x*(1 - sqrt(1-4x))/(1-4x)
F_S = x * (1 - sp.sqrt(1 - 4*x)) / (1 - 4*x)

# Expand the series up to order 16 (i.e., n = 1 to 15)
order_FS = 16
F_S_series = sp.series(F_S, x, 0, order_FS).removeO().expand()

# Extract coefficients: S(n) from generating function for n=1,2,...,15
Sackin_candidate_coeffs = [sp.nsimplify(F_S_series.coeff(x, n)) for n in range(1, order_FS)]

# -----------------------------
# Part 2: Compute Sackin Index via Recurrence
# -----------------------------
def catalan(n):
    """Return the nth Catalan number."""
    return sp.binomial(2*n, n) // (n+1)

def T_n(n):
    """
    Return the number of full binary trees with n leaves.
    T(1) = 1; for n>=2, T(n) = Catalan(n-1).
    """
    if n == 1:
        return sp.Integer(1)
    else:
        return catalan(n-1)

# Compute S(n) using the recurrence:
# S(1) = 0, and for n>=2:
# S(n) = sum_{i=1}^{n-1} [ S(i)*T(n-i) + S(n-i)*T(i) + n * T(i)*T(n-i) ]
N_max = 15
S = {1: sp.Integer(0)}  # Base case: S(1) = 0

for n in range(2, N_max+1):
    s_val = sp.Integer(0)
    for i in range(1, n):
        s_val += S[i] * T_n(n-i) + S[n-i] * T_n(i) + n * T_n(i) * T_n(n-i)
    S[n] = sp.simplify(s_val)

# Create list of S(n) values for n=1 to 15 from recurrence
Sackin_recurrence = [S[n] for n in range(1, N_max+1)]

# -----------------------------
# Part 3: Compare Generating Function Coefficients with Recurrence Values
# -----------------------------
df_comparison = pd.DataFrame({
    "n": list(range(1, N_max+1)),
    "Sackin (Generating Function)": Sackin_candidate_coeffs,
    "Sackin (Recurrence)": Sackin_recurrence,
    "Delta": [sp.nsimplify(Sackin_candidate_coeffs[i] - Sackin_recurrence[i]) for i in range(N_max)]
})
print("Comparison of Sackin Index values:")
print(df_comparison.to_string(index=False))

# -----------------------------
# Part 4: Convergence and Growth Tests
# -----------------------------
# Convert the recurrence S(n) values to float for numerical tests.
# (Skip S(1)=0 since it would cause division by zero.)
Sackin_numerical = [float(s) for s in Sackin_recurrence]
n_values = np.array(range(1, N_max+1))
# Skip the first value since S(1)=0
Sackin_nonzero = Sackin_numerical[1:]
n_nonzero = n_values[1:]

# Test 1: Ratio Test: S(n) / S(n-1) for n>=3 (using nonzero values)
ratio_test = [Sackin_nonzero[n] / Sackin_nonzero[n-1] for n in range(1, len(Sackin_nonzero))]

# Test 2: Second-Order Ratio Test: S(n) / S(n-2) for n>=4
second_order_ratio_test = [Sackin_nonzero[n] / Sackin_nonzero[n-2] for n in range(2, len(Sackin_nonzero))]

# Test 3: Root Test: S(n)^(1/n) for n>=2
root_test = [Sackin_nonzero[n-1] ** (1/n_nonzero[n-1]) for n in range(2, len(n_nonzero)+1)]

# Test 4: Difference Test: S(n) - S(n-1)
difference_test = [Sackin_nonzero[n] - Sackin_nonzero[n-1] for n in range(1, len(Sackin_nonzero))]

# Test 5: Log Growth Test: log(S(n)) - log(S(n-1))
log_growth_test = [np.log(Sackin_nonzero[n]) - np.log(Sackin_nonzero[n-1]) for n in range(1, len(Sackin_nonzero))]

# Test 6: Log-Log Growth Trend: Plot log(S(n)) vs. log(n)
log_n = np.log(n_nonzero)
log_S = np.log(Sackin_nonzero)

# -----------------------------
# Part 5: Plotting the Convergence and Growth Tests
# -----------------------------
plt.figure(figsize=(12, 10))

# Subplot 1: Ratio Test
plt.subplot(2, 3, 1)
plt.plot(n_nonzero[1:], ratio_test, 'o-', label="S(n)/S(n-1)")
plt.axhline(y=4, color='r', linestyle='--', label="Expected Limit (4)")
plt.xlabel("n")
plt.ylabel("Ratio")
plt.title("Ratio Test")
plt.legend()
plt.grid(True)

# Subplot 2: Second-Order Ratio Test
plt.subplot(2, 3, 2)
plt.plot(n_nonzero[2:], second_order_ratio_test, 's-', label="S(n)/S(n-2)")
plt.axhline(y=16, color='r', linestyle='--', label="Expected ~16 (4^2)")
plt.xlabel("n")
plt.ylabel("Second-Order Ratio")
plt.title("Second-Order Ratio Test")
plt.legend()
plt.grid(True)

# Subplot 3: Root Test
plt.subplot(2, 3, 3)
plt.plot(n_nonzero[1:], root_test, 's-', label="S(n)^(1/n)")
plt.axhline(y=4, color='r', linestyle='--', label="Expected Limit (4)")
plt.xlabel("n")
plt.ylabel("Root Value")
plt.title("Root Test")
plt.legend()
plt.grid(True)

# Subplot 4: Difference Test
plt.subplot(2, 3, 4)
plt.plot(n_nonzero[1:], difference_test, '^-', label="S(n) - S(n-1)")
plt.xlabel("n")
plt.ylabel("Difference")
plt.title("Difference Test")
plt.legend()
plt.grid(True)

# Subplot 5: Log Growth Test
plt.subplot(2, 3, 5)
plt.plot(n_nonzero[1:], log_growth_test, 'o-', label="log(S(n)) - log(S(n-1))")
plt.xlabel("n")
plt.ylabel("Log Difference")
plt.title("Log Growth Test")
plt.legend()
plt.grid(True)

# Subplot 6: Log-Log Growth Trend
plt.subplot(2, 3, 6)
plt.plot(log_n, log_S, 'bo-', label="log(S(n)) vs. log(n)")
# Fit a line to the log-log data to estimate slope (growth exponent)
slope, intercept = np.polyfit(log_n, log_S, 1)
plt.plot(log_n, slope * log_n + intercept, 'r--', label=f"Slope ~ {slope:.2f}")
plt.xlabel("log(n)")
plt.ylabel("log(S(n))")
plt.title("Log-Log Growth Trend")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# -----------------------------
# Part 6: Print Final Convergence Test Values
# -----------------------------
expanded_convergence_results = {
    "Ratio Test (Last Value, S(n)/S(n-1))": ratio_test[-1],
    "Second-Order Ratio Test (Last Value, S(n)/S(n-2))": second_order_ratio_test[-1],
    "Root Test (Last Value, S(n)^(1/n))": root_test[-1],
    "Difference Test (Last Value, S(n)-S(n-1))": difference_test[-1],
    "Log Growth Test (Last Value, log(S(n))-log(S(n-1)))": log_growth_test[-1],
    "Log-Log Growth Slope (Estimated)": slope,
    "Expected Theoretical Limit (Ratio Test)": 4
}

print("\nExpanded Convergence Test Results:")
for k, v in expanded_convergence_results.items():
    print(f"{k}: {v}")

# End of script
