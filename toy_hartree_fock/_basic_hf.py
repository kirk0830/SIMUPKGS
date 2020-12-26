# basic formulation plz see SIMUPKGS/Hartree-Fock-Roothaan.afx
# Step 1, read-in atomic coordinates

# H_{core} = T_{core} + V_{core}
# Step 2-1, integrate all kinetic parts of atomic basis T_{core}
# T_ij = <i|t|j>, i and j are basis, i or j = STO-nG

# |i> = c_1*G_1 + ...
# store |i> as a list: [(number of terms), (para pair1), (para pair2), ...]
# structure of (number of terms): (N total terms, N s-terms, N p-terms, N d-terms, ...)
# example: |i> as a 6-31G atomic basis:
# |i> = 
# [
# 10, (10 terms Gaussian function)
# 6, (6 terms s-type Gau)
# 3, (3 terms p-type Gau)
# 1, (external shell Gau)
# (C1, alpha1), (C2, alpha2), (C3, alpha3), (C4, alpha4), (C5, alpha5), (C6, alpha6),
# (C7, alpha7, 1, 0, 0), (C8, alpha8, 0, 1, 0), (C9, alpha9, 0, 0, 1), 
# (C10, alpha10, ?, ?, ?)  
# ]

# define a function makes integration between STO-nGs, support n_i != n_j

# Step 2-2, integrate all nuclear attraction parts
# V_ij = <i|V(R)|j>, i and j are basis, i or j = STO-nG

# Step 2-3, H_{core} = T_{core} + V_{core}

# Step 2-4, integrate all 4-center integrations
# ijkl-term: <ik|jl>-<ik|lj>

# Step 3, integrate all overlap integration
# S_ij = <i|j>, i and j are basis, i or j = STO-nG

# Step 4-1, diagonalize S -> s, U

# Step 4-2, sort eigenvalue of s -> s*, record number of eigenvalues to be preserved Nsp.
#           rearrange U in the same order of s* -> U*
#           and F -> F*

# Step 4-3, X = U* divides sqrt(s*) and delete last Nsp columns of X
#           then F_new = X^{-1}FX

# Step 4-4, diagonalize F_new to find C'=X^{-1}C

# Step 4-5, P=C*C

# ------Iteration start!------
# Step i-1, G_ij = sum_over_all_other_basis_by_k_and_l{P_kl * [<ik|jl>-<ik|lj>]}
#           F_ij = H_core_ij + G_ij

# Step i-2, F_new = X^{-1}FX, diagonalize F_new to find C'=X^{-1}C

# Step i-3, P=C*C
'''
C=\left( \begin{matrix}
	c_{11}&		c_{12}&		\cdots&		c_{1N}\\
	c_{21}&		c_{22}&		\cdots&		c_{2N}\\
	\vdots&		\vdots&		\ddots&		\vdots\\
	c_{N1}&		c_{N2}&		\cdots&		c_{NN}\\
\end{matrix} \right) 

P_{\mu \nu}=\sum_i{c_{\mu i}c_{i\nu}}=c_{\mu 1}c_{1\nu}+c_{\mu 2}c_{2\nu}+\cdots =CC
'''
