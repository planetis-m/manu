## Cholesky Decomposition.
##
## For a symmetric, positive definite matrix A, the Cholesky decomposition
## s an lower triangular matrix L so that A = L*L'.
##
## If the matrix is not symmetric or positive definite, the constructor
## returns a partial decomposition and sets an internal flag that may
## be queried by the isSPD() method.
import math

template newCholData() =
   newSeq(result.data, result.n)
   for i in 0 ..< result.n:
      newSeq(result.data[i], result.n)

type CholeskyDecomposition = object
   # internal array storage.
   data: seq[seq[float]]
   # matrix dimension.
   n: int
   # is symmetric and positive definite flag.
   isspd: bool

proc chol*(m: Matrix): CholeskyDecomposition =
   ## Cholesky algorithm for symmetric and positive definite matrix.
   ## Structure to access L and isspd flag.
   ## param  Square, symmetric matrix.
   result.n = m.m
   newCholData()
   result.isspd = m.n == n
   # Main loop.
   for j in 0 ..< result.n:
      var l_rowj = result.data[j]
      var d = 0.0
      for k in 0 ..< j:
         var l_rowk = result.data[k]
         var s = 0.0
         for i in 0 ..< k:
            s += l_rowk[i] * l_rowj[i]
         l_rowj[k] = s = (m.data[j][k] - s) / result.data[k][k]
         d = d + s * s
         result.isspd = result.isspd and m.data[k][j] == m.data[j][k]
      d = m.data[j][j] - d
      result.isspd = result.isspd and d > 0.0
      result.data[j][j] = sqrt(max(d, 0.0))
      for k in j + 1 ..< n:
         result.data[j][k] = 0.0

proc isSPD*(c: CholeskyDecomposition): bool =
   ## Is the matrix symmetric and positive definite?
   c.isspd

proc getL*(c: CholeskyDecomposition): Matrix =
   ## Return triangular factor.
   result.data = c.data
   result.m = c.n
   result.n = c.n

proc solve*(c: CholeskyDecomposition, b: Matrix): Matrix =
   ## Solve ``A*X = B``,
   ## returns X so that L*L'*X = B
   ## param b:  A Matrix with as many rows as A and any number of columns.
   assert(b.m == c.n, "Matrix row dimensions must agree.")
   assert(c.isspd, "Matrix is not symmetric positive definite.")

   # Copy right hand side.
   result.data = b.data
   result.m = c.n
   result.n = b.n

   # Solve L*Y = B
   for k in 0 ..< result.m:
      for j in 0 ..< result.n:
         for i in 0 ..< k:
	         result.data[k][j] -= result.data[i][j] * c.data[k][i]
	      result.data[k][j] /= c.data[k][k]

   ## Solve L'*X = Y
   for k in countdown(result.m - 1, 0):
      for j in 0 ..< result.n:
         for i in k + 1 ..< result.m:
            result.data[k][j] -= result.data[i][j] * c.data[i][k]
         result.data[k][j] /= c.data[k][k]
