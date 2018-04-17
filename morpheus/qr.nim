## QR Decomposition.
##
## For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
## orthogonal matrix Q and an n-by-n upper triangular matrix R so that
## A = Q*R.
##
## The QR decompostion always exists, even if the matrix does not have
## full rank, so the constructor will never fail.  The primary use of the
## QR decomposition is in the least squares solution of nonsquare systems
## of simultaneous linear equations.  This will fail if isFullRank()
## returns false.
import math
import "./matrix"

template newData() =
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

type QRDecomposition* = object
   # Array for internal storage of decomposition.
   data: seq[seq[float]]
   # Row and column dimensions.
   m, n: int
   # Array for internal storage of diagonal of R.
   rDiag: seq[float]

proc qr*(m: Matrix): QRDecomposition =
   ## QR Decomposition, computed by Householder reflections.
   ## Structure to access R and the Householder vectors and compute Q.
   ## parameter ``a``: Rectangular matrix

   result.data = m.data
   result.m = m.m
   result.n = m.n
   newSeq(result.rDiag, result.n)

   # Main loop.
   for k in 0 ..< result.n:
      # Compute 2-norm of k-th column without under/overflow.
      var nrm = 0.0
      for i in k ..< result.m:
         nrm = hypot(nrm, result.data[i][k])

      if nrm != 0.0:
         # Form k-th Householder vector.
         if result.data[k][k] < 0:
            nrm = -nrm
         for i in k ..< result.m:
            result.data[i][k] /= nrm
         result.data[k][k] += 1.0

         # Apply transformation to remaining columns.
         for j in k + 1 ..< result.n:
            var s = 0.0
            for i in k ..< result.m:
               s += result.data[i][k] * result.data[i][j]
            s = -s / result.data[k][k]
            for i in k ..< result.m:
               result.data[i][j] += s * result.data[i][k]
      result.rDiag[k] = -nrm

proc isFullRank*(q: QRDecomposition): bool =
   ## Is the matrix full rank?
   for j in 0 ..< q.n:
      if q.rDiag[j] == 0:
         return false
   return true

proc getH*(q: QRDecomposition): Matrix =
   ## Return the Householder vectors.
   ## return: Lower trapezoidal matrix whose columns define the reflections.
   result.m = q.m
   result.n = q.n
   newData()
   for i in 0 ..< q.m:
      for j in 0 ..< q.n:
         if i >= j:
            result.data[i][j] = q.data[i][j]

proc getR*(q: QRDecomposition): Matrix =
   ## Return the upper triangular factor.
   result.m = q.n
   result.n = q.n
   newData()
   for i in 0 ..< q.n:
      for j in 0 ..< q.n:
         if i < j:
            result.data[i][j] = q.data[i][j]
         elif i == j:
            result.data[i][j] = q.rDiag[i]
            
proc getQ*(q: QRDecomposition): Matrix =
   ## Generate and return the (economy-sized) orthogonal factor.
   result.m = q.m
   result.n = q.n
   newData()
   for k in countdown(result.n - 1, 0):
      for i in 0 ..< result.m:
         result.data[i][k] = 0.0
      result.data[k][k] = 1.0
      for j in k ..< result.n:
         if q.data[k][k] != 0:
            var s = 0.0
            for i in k ..< result.m:
               s += q.data[i][k] * result.data[i][j]
            s = -s / q.data[k][k]
            for i in k ..< result.m:
               result.data[i][j] += s * q.data[i][k]

proc solve*(q: QRDecomposition, b: Matrix): Matrix =
   ## Least squares solution of ``A*X = B``,
   ## parameter ``b``: A Matrix with as many rows as A and any number of columns.
   ## return: ``X`` that minimizes the two norm of ``Q*R*X-B``
   assert(b.m == q.m, "Matrix row dimensions must agree.")
   assert(q.isFullRank(), "Matrix is rank deficient.")

   # Copy right hand side
   result.m = q.n
   result.n = b.n
   var x = b.data

   # Compute Y = transpose(Q)*B
   for k in 0 ..< q.n:
      for j in 0 ..< b.n:
         var s = 0.0 
         for i in k ..< q.m:
            s += q.data[i][k] * x[i][j]
         s = -s / q.data[k][k]
         for i in k ..< q.m:
            x[i][j] += s * q.data[i][k]
   # Solve R*X = Y
   for k in countdown(result.n - 1, 0):
      for j in 0 ..< b.n:
         x[k][j] /= q.rDiag[k]
      for i in 0 ..< k:
         for j in 0 ..< b.n:
            x[i][j] -= x[k][j] * q.data[i][k]

   newData()
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result.data[i][j] = x[i][j]
