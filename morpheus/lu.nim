## LU Decomposition.
##
## For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
## unit lower triangular matrix L, an n-by-n upper triangular matrix U,
## and a permutation vector piv of length m so that A(piv,:) = L*U.
## If m < n, then L is m-by-m and U is m-by-n.
##
## The LU decompostion with pivoting always exists, even if the matrix is
## singular, so the constructor will never fail.  The primary use of the
## LU decomposition is in the solution of square systems of simultaneous
## linear equations.  This will fail if isNonsingular() returns false.
import "./matrix"

template newData() =
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

type LUDecomposition* = object
   # Array for internal storage of decomposition.
   lu: seq[seq[float]]
   # Row and column dimensions, and pivot sign.
   m, n, pivsign: int
   # Internal storage of pivot vector.
   piv: seq[int]

proc lu*(a: Matrix): LUDecomposition =
   ## LU Decomposition
   ## Structure to access L, U and piv.
   ## param ``a``: Rectangular matrix

   # Use a "left-looking", dot-product, Crout/Doolittle algorithm.
   result.lu = a.data
   result.m = a.m
   result.n = a.n
   newSeq(result.piv, result.m)
   for i in 0 ..< result.m:
      result.piv[i] = i
   result.pivsign = 1
   var lu_colj = newSeq[float](result.m)

   # Outer loop.
   for j in 0 ..< result.n:

      # Make a copy of the j-th column to localize references.
      for i in 0 ..< result.m:
         lu_colj[i] = result.lu[i][j]

      # Apply previous transformations.
      for i in 0 ..< result.m:
         var lu_rowi = addr result.lu[i]

         # Most of the time is spent in the following dot product.
         let kmax = min(i, j)
         var s = 0.0
         for k in 0 ..< kmax:
            s += lu_rowi[k] * lu_colj[k]
         lu_colj[i] -= s
         lu_rowi[j] = lu_colj[i]

      # Find pivot and exchange if necessary.
      var p = j
      for i in j + 1 ..< result.m:
         if abs(lu_colj[i]) > abs(lu_colj[p]):
            p = i
      if p != j:
         for k in 0 ..< result.n:
            swap(result.lu[p][k], result.lu[j][k])
         swap(result.piv[p], result.piv[j])
         result.pivsign = -result.pivsign

      # Compute multipliers.
      if j < result.m and result.lu[j][j] != 0.0:
         for i in j + 1 ..< result.m:
            result.lu[i][j] /= result.lu[j][j]

proc isNonsingular*(l: LUDecomposition): bool =
   ## Is the matrix nonsingular?
   ## returns true if U, and hence A, is nonsingular.
   for j in 0 ..< l.n:
      if l.lu[j][j] == 0.0:
         return false
   return true

proc getL*(l: LUDecomposition): Matrix =
   ## Return lower triangular factor
   result.m = l.m
   result.n = l.n
   newData()
   for i in 0 ..< l.m:
      for j in 0 ..< l.n:
         if i > j:
            result.data[i][j] = l.lu[i][j]
         elif i == j:
            result.data[i][j] = 1.0

proc getU*(l: LUDecomposition): Matrix =
   ## Return upper triangular factor
   result.m = l.n
   result.n = l.n
   newData()
   for i in 0 ..< l.n:
      for j in 0 ..< l.n:
         if i <= j:
            result.data[i][j] = l.lu[i][j]

proc getPivot*(l: LUDecomposition): seq[int] =
   ## Return pivot permutation vector
   l.piv

proc getFloatPivot*(l: LUDecomposition): seq[float] =
   ## Return pivot permutation vector as a one-dimensional double array
   result = newSeq[float](l.m)
   for i in 0 ..< l.m:
      result[i] = float(l.piv[i])

proc det*(l: LUDecomposition): float =
   ## Determinant
   assert(l.m == l.n, "Matrix must be square.")
   result = float(l.pivsign)
   for j in 0 ..< l.n:
      result *= l.lu[j][j]

proc solve*(l: LUDecomposition, b: Matrix): Matrix =
   ## Solve ``A*X = B``.
   ## parameter ``B``: A Matrix with as many rows as A and any number of columns.
   ## returns X so that ``L*U*X = B(piv,:)``
   assert(b.m == l.m, "Matrix row dimensions must agree.")
   assert(l.isNonsingular, "Matrix is singular.")

   # Copy right hand side with pivoting
   result.m = l.piv.len
   result.n = b.n
   newData()
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result.data[i][j] = b.data[l.piv[i]][j]   

   # Solve L*Y = B(piv,:)
   for k in 0 ..< l.n:
      for i in k + 1 ..< l.n:
         for j in 0 ..< result.n:
            result.data[i][j] -= result.data[k][j] * l.lu[i][k]

   # Solve U*X = Y
   for k in countdown(l.n - 1, 0):
      for j in 0 ..< result.n:
         result.data[k][j] /= l.lu[k][k]
      for i in 0 ..< k:
         for j in 0 ..< result.n:
            result.data[i][j] -= result.data[k][j] * l.lu[i][k]
