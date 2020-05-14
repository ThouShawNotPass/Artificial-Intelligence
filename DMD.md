# Dynamic Mode Decomposition

1. Convert 3D Coordinates => Hilbert => DMD
2. What is DMD?
     - X  = all atoms, all cols minus last frame
     - X' = all atoms, all cols minus first frame

3. Do SVD (Singular Valude Decomp) on X
     - Turns X into U * E * V => E is eigenvalue diagonal matrix
     - Koopman Operator
     - A * X = X' therefore
     - A * (U * E * V) = X'
