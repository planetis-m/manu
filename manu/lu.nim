## LU Decomposition
## ================
##
## For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
## unit lower triangular matrix L, an n-by-n upper triangular matrix U,
## and a permutation vector piv of length m so that ``A(piv,:) = L*U``.
## If m < n, then L is m-by-m and U is m-by-n.
##
## The LU decompostion with pivoting always exists, even if the matrix is
## singular, so the constructor will never fail. The primary use of the
## LU decomposition is in the solution of square systems of simultaneous
## linear equations. This will fail if `isNonsingular()` returns false.
import matrix

type
  LUDecomposition*[T] = object
    lu: Matrix[T] # Array for internal storage of decomposition.
    piv: seq[int] # Internal storage of pivot vector.
    pivSign: int  # Pivot sign.

proc lu*[T](a: sink Matrix[T]): LUDecomposition[T] =
  ## LU Decomposition
  ## Structure to access L, U and piv.
  ##
  ## parameter ``a``: Rectangular matrix
  # Use a "left-looking", dot-product, Crout/Doolittle algorithm.
  let m = a.m
  let n = a.n
  result.lu = a
  result.piv = newSeq[int](m)
  for i in 0 ..< m:
    result.piv[i] = i
  result.pivSign = 1
  var luColj = newSeq[float](m)
  # Outer loop.
  for j in 0 ..< n:
    # Make a copy of the j-th column to localize references.
    for i in 0 ..< m:
      luColj[i] = result.lu[i, j]
    # Apply previous transformations.
    for i in 0 ..< m:
      # Most of the time is spent in the following dot product.
      let kmax = min(i, j)
      var s = T(0)
      for k in 0 ..< kmax:
        s += result.lu[i, k] * luColj[k]
      luColj[i] -= s
      result.lu[i, j] = luColj[i]
    # Find pivot and exchange if necessary.
    var p = j
    for i in j + 1 ..< m:
      if abs(luColj[i]) > abs(luColj[p]):
        p = i
    if p != j:
      for k in 0 ..< n:
        swap(result.lu[p, k], result.lu[j, k])
      swap(result.piv[p], result.piv[j])
      result.pivSign = -result.pivSign
    # Compute multipliers.
    if j < m and result.lu[j, j] != 0:
      for i in j + 1 ..< m:
        result.lu[i, j] /= result.lu[j, j]

proc luGauss*[T](a: sink Matrix[T]): LUDecomposition[T] =
  ## LU Decomposition, computed by Gaussian elimination.
  ##
  ## This constructor computes L and U with the "daxpy"-based elimination
  ## algorithm used in LINPACK and MATLAB.
  ##
  ## - parameter ``a``: Rectangular matrix
  ## - ``return``: Structure to access L, U and piv.
  # Initialize.
  let m = a.m
  let n = a.n
  result.lu = a
  result.piv = newSeq[int](m)
  for i in 0 ..< m:
    result.piv[i] = i
  result.pivSign = 1
  # Main loop.
  for k in 0 ..< n:
    # Find pivot.
    var p = k
    for i in k + 1 ..< m:
      if abs(result.lu[i, k]) > abs(result.lu[p, k]):
        p = i
    # Exchange if necessary.
    if p != k:
      for j in 0 ..< n:
        swap(result.lu[p, j], result.lu[k, j])
      swap(result.piv[p], result.piv[k])
      result.pivSign = -result.pivSign
    # Compute multipliers and eliminate k-th column.
    if result.lu[k, k] != 0:
      for i in k + 1 ..< m:
        result.lu[i, k] /= result.lu[k, k]
        for j in k + 1 ..< n:
          result.lu[i, j] -= result.lu[i, k] * result.lu[k, j]

proc isNonsingular*[T](l: LUDecomposition[T]): bool =
  ## Is the matrix nonsingular?
  ##
  ## ``return``: true if U, and hence A, is nonsingular.
  for j in 0 ..< l.lu.n:
    if l.lu[j, j] == 0:
      return false
  return true

proc getL*[T](l: LUDecomposition[T]): Matrix[T] =
  ## Return lower triangular factor.
  result = matrixUninit[T](l.lu.m, l.lu.n)
  for i in 0 ..< l.lu.m:
    for j in 0 ..< l.lu.n:
      if i > j:
        result[i, j] = l.lu[i, j]
      elif i == j:
        result[i, j] = T(1)
      else:
        result[i, j] = T(0)

proc getU*[T](l: LUDecomposition[T]): Matrix[T] =
  ## Return upper triangular factor.
  result = matrixUninit[T](l.lu.n, l.lu.n)
  for i in 0 ..< l.lu.n:
    for j in 0 ..< l.lu.n:
      if i <= j:
        result[i, j] = l.lu[i, j]
      else:
        result[i, j] = T(0)

proc getPivot*[T](l: LUDecomposition[T]): lent seq[int] {.inline.} =
  ## Return pivot permutation vector.
  l.piv

proc getFloatPivot*[T](l: LUDecomposition[T]): seq[T] =
  ## Return pivot permutation vector as a one-dimensional double array.
  result = newSeq[T](l.lu.m)
  for i in 0 ..< l.lu.m:
    result[i] = T(l.piv[i])

proc det*[T](l: LUDecomposition[T]): T =
  ## Determinant
  assert(l.lu.m == l.lu.n, "Matrix must be square.")
  result = T(l.pivSign)
  for j in 0 ..< l.lu.n:
    result *= l.lu[j, j]

proc solve*[T](l: LUDecomposition[T], b: Matrix[T]): Matrix[T] =
  ## Solve ``A*X = B``.
  ##
  ## - parameter ``b``: A Matrix with as many rows as A and any number of columns.
  ## - ``return``: X so that ``L*U*X = B(piv,:)``
  assert(b.m == l.lu.m, "Matrix row dimensions must agree.")
  assert(l.isNonsingular, "Matrix is singular.")
  # Copy right hand side with pivoting
  result = b[l.piv, 0 ..< b.n]
  # Solve L*Y = B(piv,:)
  for k in 0 ..< l.lu.n:
    for i in k + 1 ..< l.lu.n:
      for j in 0 ..< b.n:
        result[i, j] -= result[k, j] * l.lu[i, k]
  # Solve U*X = Y
  for k in countdown(l.lu.n - 1, 0):
    for j in 0 ..< b.n:
      result[k, j] /= l.lu[k, k]
    for i in 0 ..< k:
      for j in 0 ..< b.n:
        result[i, j] -= result[k, j] * l.lu[i, k]
