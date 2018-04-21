## Cholesky Decomposition.
##
## For a symmetric, positive definite matrix A, the Cholesky decomposition
## is an lower triangular matrix L so that A = L*L'.
##
## If the matrix is not symmetric or positive definite, the constructor
## returns a partial decomposition and sets an internal flag that may
## be queried by the isSPD() method.
import math
import "./matrix"

type CholeskyDecomposition* = object
   # triangular factor
   l: Matrix
   # is symmetric and positive definite flag.
   isspd: bool

proc chol*(a: Matrix): CholeskyDecomposition =
   ## Cholesky algorithm for symmetric and positive definite matrix.
   ## Structure to access L and isspd flag.
   ## parameter ``m``: Square, symmetric matrix.
   let n = a.m
   result.l = matrix(n, n)
   result.isspd = a.n == a.m
   # Main loop.
   for j in 0 ..< n:
      var lRowj = addr result.l.mgetRow(j)
      var d = 0.0
      for k in 0 ..< j:
         var lRowk = addr result.l.mgetRow(k)
         var s = 0.0
         for i in 0 ..< k:
            s += lRowk[i] * lRowj[i]
         s = (a[j, k] - s) / result.l[k, k]
         lRowj[k] = s
         d = d + s * s
         result.isspd = result.isspd and a[k, j] == a[j, k]
      d = a[j, j] - d
      result.isspd = result.isspd and d > 0.0
      result.l[j, j] = sqrt(max(d, 0.0))
#       for k in j + 1 ..< n:
#          result.l[j, k] = 0.0

proc isSPD*(c: CholeskyDecomposition): bool {.inline.} =
   ## Is the matrix symmetric and positive definite?
   c.isspd

proc getL*(c: CholeskyDecomposition): Matrix {.inline.} =
   ## Return triangular factor.
   c.l

proc solve*(c: CholeskyDecomposition, b: Matrix): Matrix =
   ## Solve ``A*X = B``,
   ## parameter ``b``:  A Matrix with as many rows as A and any number of columns.
   ## return: X so that ``L*L'*X = B``
   let n = c.l.m
   let nx = b.n
   assert(b.m == n, "Matrix row dimensions must agree.")
   assert(c.isspd, "Matrix is not symmetric positive definite.")
   # Copy right hand side.
   result = b
   # Solve L*Y = B
   for k in 0 ..< n:
      for j in 0 ..< nx:
         for i in 0 ..< k:
            result[k, j] -= result[i, j] * c.l[k, i]
         result[k, j] /= c.l[k, k]
   # Solve L'*X = Y
   for k in countdown(n - 1, 0):
      for j in 0 ..< nx:
         for i in k + 1 ..< n:
            result[k, j] -= result[i, j] * c.l[i, k]
         result[k, j] /= c.l[k, k]
