## Morpheus - Nim Matrix library
## =============================
##
## The morpheus module provides the fundamental operations of numerical
## linear algebra. Various constructors create Matrices from two dimensional
## arrays of double precision floating point numbers. Various "gets" and
## "sets" provide access to submatrices and matrix elements.  Several methods 
## implement basic matrix arithmetic, including matrix addition and
## multiplication, matrix norms, and element-by-element array operations.
## Methods for reading and printing matrices are also included. All the
## operations in this version of the Matrix object involve real matrices.
## Complex matrices may be handled in a future version.
##
## Five fundamental matrix decompositions, which consist of pairs or triples
## of matrices, permutation vectors, and the like, produce results in five
## decomposition classes. These decompositions are accessed by the morpheus
## module to compute solutions of simultaneous linear equations, determinants,
## inverses and other matrix functions. The five decompositions are:
##
## - Cholesky Decomposition of symmetric, positive definite matrices.
## - LU Decomposition of rectangular matrices.
## - QR Decomposition of rectangular matrices.
## - Singular Value Decomposition of rectangular matrices.
## - Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
##
## **Example of use:**
##
## .. code-block:: nim
##    import morpheus
##    # Solve a linear system A x = b and compute the residual norm, ||b - A x||.
##    let vals = @[@[1.0, 2.0, 3.0], @[4.0, 5.0, 6.0], @[7.0, 8.0, 10.0]]
##    let A = matrix(vals)
##    let b = randMatrix(3, 1)
##    let x = A.solve(b)
##    let r = A * x - b
##    let rnorm = r.normInf()
import "./morpheus" / [matrix, cholesky, qr, lu, svd, eigen]
export matrix, cholesky, qr, lu, svd, eigen

proc norm2*(m: Matrix): float =
   ## Two norm
   ##
   ## ``return``: maximum singular value.
   svd(m).norm2()

proc solve*(a, b: Matrix): Matrix =
   ## Solve ``A*X = B``,
   ##
   ## - parameter ``b``: the right hand side
   ## - ``return``: solution if A is square, least squares solution otherwise.
   if a.m == a.n:
      lu(a).solve(b)
   else:
      qr(a).solve(b)

proc solveTranspose*(a, b: Matrix): Matrix =
   ## Solve ``X*A = B``, which is also ``A'*X' = B'``,
   ##
   ## - parameter ``b``: the right hand side
   ## - ``return``: solution if A is square, least squares solution otherwise.
   transpose(a).solve(b.transpose())

proc inverse*(m: Matrix): Matrix =
   ## Matrix inverse or pseudoinverse,
   ##
   ## ``return``: inverse(A) if A is square, pseudoinverse otherwise.
   solve(m, identity(m.m, m.m))

proc det*(m: Matrix): float =
   ## Matrix determinant
   lu(m).det()

proc rank*(m: Matrix): int =
   ## Matrix rank,
   ##
   ## ``return``: effective numerical rank, obtained from SVD.
   svd(m).rank()

proc cond*(m: Matrix): float =
   ## Matrix condition (2 norm),
   ##
   ## ``return``: ratio of largest to smallest singular value.
   svd(m).cond()
