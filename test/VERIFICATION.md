# Verification Methodology for linear-massiv

This document describes the numerical verification methodology used in the
linear-massiv test suite. The approach follows established best practices
from the LAPACK project, Higham's error analysis framework, and the
numerical linear algebra literature.

## 1. Backward Error Analysis

The foundation of our verification is *backward error analysis*, as
pioneered by Wilkinson and systematized by Higham (2002). Rather than
measuring the forward error (distance from the true solution), we measure
the *backward error*: how much must the input be perturbed to make the
computed output exact?

For a computed solution x-hat to Ax = b, the normwise backward error is:

    eta(x-hat) = min { eps : (A + dA) x-hat = b + db,
                       ||dA|| <= eps ||A||, ||db|| <= eps ||b|| }

This is computable via the Rigaud formula (Higham, 2002, Theorem 7.1):

    eta(x-hat) = ||r|| / (||A|| ||x-hat|| + ||b||)

where r = b - A x-hat is the residual.

An algorithm is *backward stable* if it produces eta = O(n * eps), where
eps is machine epsilon and n is the problem dimension.

## 2. Machine Epsilon

IEEE 754 double precision unit roundoff:

    eps = 2^(-52) = 2.220446049250313e-16

This is the smallest value such that fl(1 + eps) > 1 in double precision
floating-point arithmetic. All error bounds in this test suite are
expressed in terms of eps.

## 3. Scaled Residual Formulas

The test suite uses the following scaled residuals, matching the LAPACK
testing methodology (LAWN 41; Anderson et al., 1999):

### Linear System Solve (Ax = b)

    ||A x-hat - b||_inf / (||A||_inf ||x-hat||_inf + ||b||_inf)

Per Higham (2002), Section 7.1. A backward-stable solver achieves
this quantity at O(n * eps).

### Eigenvalue Problem (Av = lambda v)

    ||Av - lambda v||_2 / (||A||_F ||v||_2)

Per Higham (2002), Section 14.1.

### QR Factorization (A = QR)

    ||A - QR||_F / (||A||_F * n * eps)

Per LAPACK testing (LAWN 41). A value less than O(1) (typically
< 10-100) indicates backward stability.

### Orthogonality (Q^T Q = I)

    ||Q^T Q - I||_F / (n * eps)

Per LAPACK testing. Values < O(1) mean orthogonality is preserved
to within rounding error.

### SVD (A = U Sigma V^T)

    ||A - U Sigma V^T||_F / (||A||_F * max(m,n) * eps)

### Cholesky (A = G G^T)

    ||A - G G^T||_F / (||A||_F * n * eps)

### LU with Partial Pivoting (PA = LU)

    ||PA - LU||_F / (||A||_F * n * eps)

## 4. Condition-Number-Aware Tolerances

For a linear system Ax = b, the forward error satisfies (Higham, 2002,
Theorem 7.2):

    ||x-hat - x|| / ||x|| <= kappa(A) * eta

where kappa(A) = ||A|| ||A^{-1}|| is the condition number and eta is
the backward error. A backward-stable solver achieves eta ~ n * eps,
so the expected forward error is:

    O(kappa * n * eps)

For test matrices with known condition numbers (e.g., from the randsvd
construction), we set tolerances proportionally to kappa * n * eps.
For general random matrices, we estimate kappa via SVD.

## 5. Test Matrix Catalog

### Hilbert Matrix

    H(i,j) = 1 / (i + j + 1)

- Symmetric positive definite
- Condition number grows exponentially: kappa_2(H_5) ~ 4.8e5,
  kappa_2(H_10) ~ 1.6e13
- Inverse has integer entries
- Tests stability of solvers under extreme ill-conditioning

Reference: Higham (2002), Section 28.1; GVL4, p. 128.

### Wilkinson Matrix

    W(i,i) = |i - floor(n/2)|,  W(i,i+1) = W(i+1,i) = 1

- Symmetric tridiagonal
- Has near-degenerate eigenvalue pairs
- Tests resolution of close eigenvalues by iterative solvers

Reference: Higham (2002), Section 28.6; Wilkinson (1965).

### Frank Matrix

    F(i,j) = n - max(i,j)  for j >= i-1,  0 otherwise

- Upper Hessenberg
- det(F) = 1
- Known real positive eigenvalues
- Ill-conditioned eigenvalues

Reference: Higham (2002), Section 28.5; Frank (1958).

### Hadamard Matrix (Sylvester/Walsh Construction)

    H(i,j) = (-1)^popcount(i AND j)

- Entries are all +/-1
- For power-of-2 sizes: H^T H = n I (orthogonal up to scaling)
- Tests preservation of orthogonality

Reference: Higham (2002), Section 28.3.

### Randsvd Construction

    A = U diag(sigma_1, ..., sigma_n) V^T

where U, V are random orthogonal matrices obtained from QR factorization
of random matrices, and the singular values are geometrically spaced:

    sigma_i = kappa^(-i/(n-1))

giving sigma_1 = 1 and sigma_n = 1/kappa, so cond_2(A) = kappa exactly.

Reference: Higham (2002), Section 28.3; Fasi & Higham (2021).

### Clustered Eigenvalue Matrices

    A = Q diag(c + i * 10^{-6}) Q^T

where Q is a random orthogonal matrix and c is the cluster center.
All eigenvalues are within a tiny interval, stressing iterative
eigenvalue solvers that must resolve near-degenerate pairs.

## 6. Property-Based vs Known-Answer Testing

The test suite uses both approaches:

**Known-answer tests** verify computed results against analytically
known solutions (e.g., eigenvalues of diagonal matrices, det of
identity). These provide definitive correctness checks for specific
cases.

**Property-based tests** (via QuickCheck) verify mathematical
invariants that must hold for all inputs:
- Reconstruction: A = QR, PA = LU, A = GG^T, A = U Sigma V^T
- Orthogonality: Q^T Q = I
- Algebraic identities: det(A) = product(eigenvalues), (AB)C = A(BC)
- Norm properties: triangle inequality, submultiplicativity

The combination provides both depth (specific known cases) and
breadth (random exploration of the input space).

## 7. References

1. Higham, N. J. (2002). *Accuracy and Stability of Numerical
   Algorithms*, 2nd ed., SIAM. ISBN 0-89871-521-0.
   - Chapter 1: Principles of finite precision computation
   - Chapter 7: Norms and error analysis for linear systems
   - Chapter 14: Eigenvalue problems
   - Chapter 28: Test Matrix Toolbox

2. Golub, G. H. & Van Loan, C. F. (2013). *Matrix Computations*,
   4th ed., Johns Hopkins University Press.
   - Chapter 2: Matrix analysis (norms, perturbation theory)
   - Section 2.7: Finite precision matrix computations

3. Anderson, E. et al. (1999). *LAPACK Users' Guide*, 3rd ed., SIAM.
   - Chapter 10: Testing and timing

4. Demmel, J. W. (1997). *Applied Numerical Linear Algebra*, SIAM.
   ISBN 0-89871-389-7.

5. LAPACK Working Note 9: Demmel, J. W. & Dongarra, J. (1989).
   "A Test Matrix Generation Suite." University of Tennessee.
   URL: https://www.netlib.org/lapack/lawnspdf/lawn09.pdf

6. LAPACK Working Note 41: Anderson, E. et al. (1999).
   "Installation Guide for LAPACK."
   URL: https://icl.utk.edu/~mgates3/docs/lawn41.pdf

7. Higham, N. J. (1997). "Testing Linear Algebra Software," in
   *Quality of Numerical Software*, Chapman and Hall, pp. 109-122.
   URL: https://nhigham.com/wp-content/uploads/2023/10/high97t.pdf

8. Fasi, M. & Higham, N. J. (2021). "Generating Extreme-Scale
   Matrices With Specified Singular Values or Condition Number,"
   SIAM J. Sci. Comput. 43(5).
   URL: https://epubs.siam.org/doi/10.1137/20M1327938

9. Kellner, T., Higham, N. J. & Mikaitis, M. (2022). "Generation of
   test matrices with specified eigenvalues using floating-point
   arithmetic," Numerical Algorithms.
   URL: https://link.springer.com/article/10.1007/s11075-021-01186-7

10. Appel, A. W. et al. (2023). "LAProof: A Library of Formal Proofs
    of Accuracy and Correctness for Linear Algebra Programs,"
    IEEE ARITH 2023.
    URL: https://www.cs.princeton.edu/~appel/papers/LAProof.pdf

11. Claessen, K. & Hughes, J. (2000). "QuickCheck: A Lightweight Tool
    for Random Testing of Haskell Programs," ICFP 2000.

12. Davis, T. A. & Hu, Y. (2011). "The University of Florida Sparse
    Matrix Collection," ACM TOMS 38(1).

13. Duff, I. S., Grimes, R. G. & Lewis, J. G. (1989). "Sparse matrix
    test problems," ACM TOMS 15(1).
