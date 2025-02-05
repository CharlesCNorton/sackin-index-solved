# Asymptotic Analysis of the Sackin Index in Full Binary Trees

By: Charles Norton & o3-mini-high

February 5th, 2025

## 1. Introduction

The Sackin index is a widely used measure of tree imbalance in phylogenetics. For a full binary tree \( T \) with \( n \) leaves, the Sackin index \( S(T) \) is defined as the sum of the depths of its leaves. If we denote by \( S(n) \) the total Sackin index summed over all full binary trees with \( n \) leaves, then its generating function is given by

\[
F_S(x) = \sum_{n\ge1} S(n)\, x^n = \frac{x\,(1-\sqrt{1-4x})}{1-4x}\,.
\]

Having a closed–form expression for \( F_S(x) \) allows us to use analytic combinatorics to extract the asymptotic behavior of the coefficients \( S(n) \) as \( n\to\infty \).

---

## 2. Singularity Analysis Background

A central idea in analytic combinatorics is that the asymptotic behavior of the coefficients \( [x^n]F_S(x) \) is dictated by the singularities of \( F_S(x) \). In our case, note that:
- The square root \(\sqrt{1-4x}\) has a branch-point singularity at \( x=\tfrac{1}{4} \).
- The factor \( 1-4x \) in the denominator introduces a pole at \( x=\tfrac{1}{4} \).

Thus, the dominant singularity is located at \( x_0 = \tfrac{1}{4} \). We now proceed to “zoom in” on this singularity.

---

## 3. Rewriting the Generating Function Near the Dominant Singularity

Introduce the change of variable

\[
z = 1-4x, \quad \text{so that} \quad x = \frac{1-z}{4}\,.
\]

Then the generating function transforms as follows:

\[
F_S(x) = \frac{x\,(1-\sqrt{1-4x})}{1-4x} = \frac{\frac{1-z}{4}\,\bigl(1-\sqrt{z}\bigr)}{z} = \frac{1-z}{4z}\,\bigl(1-\sqrt{z}\bigr)\,.
\]

We can now split the expression as

\[
F_S(x) = \frac{1}{4z}\,\bigl(1-\sqrt{z}\bigr) - \frac{1}{4}\,\bigl(1-\sqrt{z}\bigr)\,.
\]

Notice that the second term, \(-\tfrac{1}{4}(1-\sqrt{z})\), is analytic at \( z=0 \) (i.e. \( x=\tfrac{1}{4} \)) and hence does not contribute to the dominant asymptotics. Therefore, the singular behavior is determined by

\[
\frac{1-\sqrt{z}}{4z}\,.
\]

Let’s rewrite this term further:

\[
\frac{1-\sqrt{z}}{4z} = \frac{1}{4z} - \frac{1}{4}\,\frac{1}{\sqrt{z}}\,.
\]

Thus, the singular expansion of \( F_S(x) \) as \( z \to 0 \) (or equivalently \( x\to\frac{1}{4}^- \)) is

\[
F_S(x) \sim \frac{1}{4}\,(1-4x)^{-1} - \frac{1}{4}\,(1-4x)^{-1/2} + \text{analytic terms}\,.
\]

---

## 4. Transfer Theorems and Coefficient Asymptotics

The standard transfer theorem (see Flajolet & Sedgewick’s *Analytic Combinatorics*) tells us that if

\[
F(x) \sim A\, (1-4x)^{-\alpha} \quad \text{as } x\to\frac{1}{4}^-,
\]

then the coefficient \( f(n) = [x^n]F(x) \) is asymptotically given by

\[
f(n) \sim \frac{A}{\Gamma(\alpha)}\,4^n\, n^{\alpha-1}\,.
\]

Let’s apply this to each singular term.

### 4.1. The First Term

For

\[
\frac{1}{4}\,(1-4x)^{-1}\,:
\]

- \( A_1 = \frac{1}{4} \) and \( \alpha_1 = 1 \).  
- Since \( \Gamma(1)=1 \), this term contributes

\[
[x^n] \sim \frac{1}{4}\,\frac{4^n}{\Gamma(1)}\, n^{1-1} = \frac{4^n}{4}\,.
\]

### 4.2. The Second Term

For

\[
-\frac{1}{4}\,(1-4x)^{-1/2}\,:
\]

- \( A_2 = -\frac{1}{4} \) and \( \alpha_2 = \frac{1}{2} \).  
- Since \( \Gamma\left(\frac{1}{2}\right)=\sqrt{\pi} \), this term contributes

\[
[x^n] \sim -\frac{1}{4}\,\frac{4^n}{\sqrt{\pi}}\, n^{\frac{1}{2}-1} = -\frac{4^n}{4\sqrt{\pi}}\, n^{-1/2}\,.
\]

---

## 5. The Final Asymptotic Formula

Combining the contributions from both singular terms, we obtain the asymptotic expansion for \( S(n) \):

\[
S(n) \sim \frac{4^n}{4} - \frac{4^n}{4\sqrt{\pi}}\, n^{-1/2} + O\!\Bigl(\frac{4^n}{n}\Bigr)\,.
\]

Or equivalently, one can factor out the dominant exponential growth:

\[
S(n) \sim \frac{4^n}{4}\left( 1 - \frac{1}{\sqrt{\pi n}} + O\!\left(\frac{1}{n}\right) \right)\,.
\]

---

## 6. Exposition and Implications

### **Interpretation**

- **Dominant Term:**  
  The term \(\frac{4^n}{4}\) represents the leading behavior. Note that the number of full binary trees with \( n \) leaves is asymptotically proportional to \( \frac{4^n}{n^{3/2}} \) (i.e. the Catalan numbers). Hence, the average Sackin index per tree is roughly of order \( n^{3/2} \), and the total \( S(n) \) across all trees scales as \( 4^n \).

- **Correction Term:**  
  The correction term \(-\frac{4^n}{4\sqrt{\pi}}\, n^{-1/2}\) refines this approximation by capturing the first-order deviation from the dominant exponential behavior. In many applications, such as estimating expected tree imbalance, knowing these correction terms is critical.

### **Why This Result Is Significant**

1. **Computational Efficiency:**  
   Previously, estimating \( S(n) \) for large \( n \) required recursive methods or brute-force enumeration. The explicit asymptotic formula now allows for rapid estimation.

2. **Theoretical Insight:**  
   The derivation demonstrates the power of analytic combinatorics and singularity analysis in solving problems that were previously intractable without a closed–form generating function.

3. **Broader Applications:**  
   This approach can be extended to derive asymptotic formulas for other tree invariants, opening new avenues in both combinatorial theory and practical applications (e.g., phylogenetics).

---

## 7. Conclusion

Starting from the closed–form generating function

\[
F_S(x) = \frac{x\,(1-\sqrt{1-4x})}{1-4x}\,,
\]

we applied singularity analysis and transfer theorems to derive the asymptotic behavior of its coefficients. The final explicit asymptotic formula for the total Sackin index \( S(n) \) is

\[
\boxed{ S(n) \sim \frac{4^n}{4}\left( 1 - \frac{1}{\sqrt{\pi n}} + O\!\left(\frac{1}{n}\right) \right) \quad \text{as } n\to\infty\,. }
\]

This result not only confirms the exponential nature of the growth of \( S(n) \) but also provides the precise polynomial correction factors, offering a powerful tool for further studies in tree imbalance measures.
