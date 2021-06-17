## Cholesky Decomposition
## ======================
##
## For a symmetric, positive definite matrix A, the Cholesky decomposition
## is an lower triangular matrix L so that ``A = L*L'``.
##
## If the matrix is not symmetric or positive definite, the constructor
## returns a partial decomposition and sets an internal flag that may
## be queried by the `isSpd()` method.
import matrix, std/math

type
  CholeskyDecomposition*[T] = object
    l: Matrix[T] # triangular factor
    isSpd: bool  # is symmetric and positive definite flag.

proc chol*[T](a: Matrix[T]): CholeskyDecomposition[T] =
  ## Cholesky algorithm for symmetric and positive definite matrix.
  ## Structure to access L and isSpd flag.
  ##
  ## parameter ``m``: Square, symmetric matrix.
  # assert(a.n == a.m and a.n >= 1)
  let n = a.m
  result.l = matrixUninit[T](n, n)
  result.isSpd = a.n == a.m
  # Main loop.
  for j in 0 ..< n:
    var d = T(0.0)
    for k in 0 ..< j:
      var s = T(0.0)
      for i in 0 ..< k:
        s += result.l[k, i] * result.l[j, i]
      s = (a[j, k] - s) / result.l[k, k]
      result.l[j, k] = s
      d = d + s * s
      result.isSpd = result.isSpd and a[k, j] == a[j, k]
    d = a[j, j] - d
    result.isSpd = result.isSpd and d > 0.0
    result.l[j, j] = sqrt(max(d, T(0.0)))
    for k in j + 1 ..< n:
      result.l[j, k] = T(0.0)

proc isSpd*[T](c: CholeskyDecomposition[T]): bool {.inline.} =
  ## Is the matrix symmetric and positive definite?
  c.isSpd

proc getL*[T](c: CholeskyDecomposition[T]): lent Matrix[T] {.inline.} =
  ## Return triangular factor.
  c.l

proc solve*[T](c: CholeskyDecomposition[T], b: sink Matrix[T]): Matrix[T] =
  ## Solve ``A*X = B``
  ##
  ## - parameter ``b``:  A Matrix with as many rows as A and any number of columns.
  ## - ``return``: X so that ``L*L'*X = B``
  assert(b.m == c.l.n, "Matrix row dimensions must agree.")
  assert(c.isSpd, "Matrix is not symmetric positive definite.")
  # Copy right hand side.
  result = b
  # Solve L*Y = B
  for k in 0 ..< c.l.m:
    for j in 0 ..< result.n:
      for i in 0 ..< k:
        result[k, j] -= result[i, j] * c.l[k, i]
      result[k, j] /= c.l[k, k]
  # Solve L'*X = Y
  for k in countdown(c.l.m - 1, 0):
    for j in 0 ..< result.n:
      for i in k + 1 ..< c.l.m:
        result[k, j] -= result[i, j] * c.l[i, k]
      result[k, j] /= c.l[k, k]
