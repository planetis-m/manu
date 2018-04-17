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

template newData() =
   newSeq(result.data, result.n)
   for i in 0 ..< result.n:
      newSeq(result.data[i], result.n)

type CholeskyDecomposition* = object
   # Array for internal storage of decomposition.
   data: seq[seq[float]]
   # Row and column dimension (square matrix).
   n: int
   # Symmetric and positive definite flag.
   isspd: bool

proc chol*(m: Matrix): CholeskyDecomposition =
   ## Cholesky algorithm for symmetric and positive definite matrix.
   ## Structure to access L and isspd flag.
   ## parameter ``m``: Square, symmetric matrix.
   result.n = m.m
   newData()
   result.isspd = m.n == m.m
   # Main loop.
   for j in 0 ..< result.n:
      var lRowj = addr result.data[j]
      var d = 0.0
      for k in 0 ..< j:
         var lRowk = addr result.data[k]
         var s = 0.0
         for i in 0 ..< k:
            s += lRowk[i] * lRowj[i]
         s = (m.data[j][k] - s) / result.data[k][k]
         lRowj[k] = s
         d = d + s * s
         result.isspd = result.isspd and m.data[k][j] == m.data[j][k]
      d = m.data[j][j] - d
      result.isspd = result.isspd and d > 0.0
      result.data[j][j] = sqrt(max(d, 0.0))
      # for k in j + 1 ..< result.n:
      #    result.data[j][k] = 0.0

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
   ## parameter ``b``:  A Matrix with as many rows as A and any number of columns.
   ## return: X so that ``L*L'*X = B``
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

   # Solve L'*X = Y
   for k in countdown(result.m - 1, 0):
      for j in 0 ..< result.n:
         for i in k + 1 ..< result.m:
            result.data[k][j] -= result.data[i][j] * c.data[i][k]
         result.data[k][j] /= c.data[k][k]
