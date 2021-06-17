## Manu - Nim Matrix library
## =========================
##
## The manu module provides the fundamental operations of numerical
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
## decomposition classes. These decompositions are accessed by the manu
## module to compute solutions of simultaneous linear equations, determinants,
## inverses and other matrix functions. The five decompositions are:
##
## - Cholesky Decomposition of symmetric, positive definite matrices.
## - LU Decomposition of rectangular matrices.
## - QR Decomposition of rectangular matrices.
## - Singular Value Decomposition of rectangular matrices.
## - Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
runnableExamples:
  # Solve a linear system A x = b and compute the residual norm, ||b - A x||.
  let vals = @[@[1.0, 2, 3], @[4.0, 5, 6], @[7.0, 8, 10]]
  let A = matrix(vals)
  let b = randMatrix64(3, 1)
  let x = A.solve(b)
  let r = A * x - b
  let rnorm = r.normInf()
  echo("x =\n", x)
  echo("residual norm = ", rnorm)

import manu / [matrix, cholesky, qr, lu, svd, eigen]
export matrix, cholesky, qr, lu, svd, eigen

proc norm2*[T](m: sink Matrix[T]): T =
  ## Two norm
  ##
  ## ``return``: maximum singular value.
  svd(m).norm2()

proc solve*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
  ## Solve ``A*X = B``
  ##
  ## - parameter ``b``: the right hand side
  ## - ``return``: solution if A is square, least squares solution otherwise.
  if a.m == a.n:
    lu(a).solve(b)
  else:
    qr(a).solve(b)

proc solveTranspose*[T](a, b: Matrix[T]): Matrix[T] =
  ## Solve ``X*A = B``, which is also ``A'*X' = B'``
  ##
  ## - parameter ``b``: the right hand side
  ## - ``return``: solution if A is square, least squares solution otherwise.
  transpose(a).solve(b.transpose())

proc inverse*[T](m: sink Matrix[T]): Matrix[T] =
  ## Matrix inverse or pseudoinverse
  ##
  ## ``return``: inverse(A) if A is square, pseudoinverse otherwise.
  let id = identity[T](m.m)
  solve(m, id)

proc det*[T](m: sink Matrix[T]): T =
  ## Matrix determinant
  lu(m).det()

proc rank*[T](m: sink Matrix[T]): int =
  ## Matrix rank
  ##
  ## ``return``: effective numerical rank, obtained from SVD.
  svd(m).rank()

proc cond*[T](m: sink Matrix[T]): T =
  ## Matrix condition (2 norm)
  ##
  ## ``return``: ratio of largest to smallest singular value.
  svd(m).cond()
