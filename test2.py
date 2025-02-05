"""
Rigorous Unit Tests for the Closed–Form Generating Function of the Sackin Index

By: Charles Norton and o3-mini-high

The Sackin index of a full binary tree is defined as the sum of the depths of its leaves (with the root at depth 0).
The candidate generating function is:
    F_S(x) = x*(1 - sqrt(1-4*x))/(1-4*x)

We verify that:
  1. The series expansion coefficients for F_S(x) (for n = 1..15)
     match the values computed using the recurrence:
         S(1) = 0,
         S(n) = sum_{i=1}^{n-1}[ S(i)*T(n-i) + S(n-i)*T(i) + n*T(i)*T(n-i) ]
  2. For small n (n = 1..8), a brute-force tree generation method gives the same Sackin index totals.
"""

import unittest
import sympy as sp

class TestSackinIndex(unittest.TestCase):
    # --- Helper Methods ---
    def catalan(self, n):
        """Return the nth Catalan number.
           For full binary trees, T(1)=1 and for n>=2, T(n)=Catalan(n-1).
        """
        if n == 1:
            return sp.Integer(1)
        else:
            return sp.binomial(2*(n-1), n-1) // n

    def compute_T(self, n):
        """Number of full binary trees with n leaves."""
        return self.catalan(n)

    def compute_sackin_recurrence(self, N_max):
        """
        Computes S(n) using the recurrence:
          S(1)=0,
          S(n)=∑[ S(i)*T(n-i) + S(n-i)*T(i) + n*T(i)*T(n-i) ], for i=1 to n-1.
        Returns a dictionary mapping n to S(n) for n=1..N_max.
        """
        S = {1: sp.Integer(0)}
        T_cache = {n: sp.Integer(self.compute_T(n)) for n in range(1, N_max+1)}
        for n in range(2, N_max+1):
            total = sp.Integer(0)
            for i in range(1, n):
                total += S[i] * T_cache[n-i] + S[n-i] * T_cache[i] + n * T_cache[i] * T_cache[n-i]
            S[n] = sp.simplify(total)
        return S

    def candidate_coefficients(self, N_max):
        """
        Expands F_S(x) = x*(1 - sqrt(1-4*x))/(1-4*x) as a power series up to order N_max+1
        and returns a list of coefficients for x^n, n=1..N_max.
        """
        x = sp.symbols('x')
        F_S_candidate = x * (1 - sp.sqrt(1-4*x)) / (1-4*x)
        series_expansion = sp.series(F_S_candidate, x, 0, N_max+1).removeO().expand()
        coeffs = [sp.nsimplify(series_expansion.coeff(x, n)) for n in range(1, N_max+1)]
        return coeffs

    def generate_full_binary_trees(self, n):
        """
        Recursively generates all full binary trees with n leaves.
        A leaf is represented as the string "L", and an internal node as a tuple (left, right).
        """
        if n == 1:
            return ["L"]
        trees = []
        for i in range(1, n):
            left_subtrees = self.generate_full_binary_trees(i)
            right_subtrees = self.generate_full_binary_trees(n - i)
            for left in left_subtrees:
                for right in right_subtrees:
                    trees.append((left, right))
        return trees

    def sackin_index(self, tree, depth=0):
        """
        Computes the Sackin index of a given tree.
        For a leaf ("L"), returns the current depth.
        For an internal node (tuple), returns the sum of the Sackin indices of its subtrees.
        """
        if tree == "L":
            return depth
        else:
            left, right = tree
            return self.sackin_index(left, depth + 1) + self.sackin_index(right, depth + 1)

    def brute_force_sackin(self, n):
        """
        For a given n, generates all full binary trees with n leaves,
        computes each tree's Sackin index, and returns the total sum S(n).
        """
        trees = self.generate_full_binary_trees(n)
        total = sp.Integer(0)
        for tree in trees:
            total += self.sackin_index(tree)
        return total

    # --- Test Cases ---
    def test_candidate_vs_recurrence(self):
        """Verify that the candidate generating function's coefficients match the recurrence values."""
        N_max = 15
        rec_S = self.compute_sackin_recurrence(N_max)
        cand_coeffs = self.candidate_coefficients(N_max)
        for n in range(1, N_max+1):
            with self.subTest(n=n):
                self.assertEqual(rec_S[n], cand_coeffs[n-1],
                                 f"Mismatch at n={n}: recurrence {rec_S[n]}, candidate {cand_coeffs[n-1]}")

    def test_candidate_vs_bruteforce(self):
        """For small n, verify that brute-force tree generation yields the same Sackin index totals as the candidate series."""
        max_n = 8
        cand_coeffs = self.candidate_coefficients(max_n)
        for n in range(1, max_n+1):
            with self.subTest(n=n):
                brute_val = self.brute_force_sackin(n)
                self.assertEqual(cand_coeffs[n-1], brute_val,
                                 f"Mismatch at n={n}: candidate {cand_coeffs[n-1]}, brute force {brute_val}")

if __name__ == '__main__':
    unittest.main()
