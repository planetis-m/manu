## QR Decomposition.
##
## For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
## orthogonal matrix Q and an n-by-n upper triangular matrix R so that
## A = Q*R.
##
## The QR decompostion always exists, even if the matrix does not have
## full rank, so the constructor will never fail. The primary use of the
## QR decomposition is in the least squares solution of nonsquare systems
## of simultaneous linear equations. This will fail if isFullRank()
## returns false.
import "./matrix", math

type
   QRDecomposition* = object
      qr: Matrix # Array for internal storage of decomposition.
      rDiag: seq[float] # Array for internal storage of diagonal of R.

proc qr*(a: sink Matrix): QRDecomposition =
   ## QR Decomposition, computed by Householder reflections.
   ## Structure to access R and the Householder vectors and compute Q.
   ##
   ## parameter ``a``: Rectangular matrix
   let m = a.m
   let n = a.n
   result.qr = a
   result.rDiag = newSeq[float](n)
   # Main loop.
   for k in 0 ..< n:
      # Compute 2-norm of k-th column without under/overflow.
      var nrm = 0.0
      for i in k ..< m:
         nrm = hypot(nrm, result.qr[i, k])
      if nrm != 0.0:
         # Form k-th Householder vector.
         if result.qr[k, k] < 0:
            nrm = -nrm
         for i in k ..< m:
            result.qr[i, k] /= nrm
         result.qr[k, k] += 1.0
         # Apply transformation to remaining columns.
         for j in k + 1 ..< n:
            var s = 0.0
            for i in k ..< m:
               s += result.qr[i, k] * result.qr[i, j]
            s = -s / result.qr[k, k]
            for i in k ..< m:
               result.qr[i, j] += s * result.qr[i, k]
      result.rDiag[k] = -nrm

proc isFullRank*(q: QRDecomposition): bool =
   ## Is the matrix full rank?
   for d in q.rDiag:
      if d == 0.0:
         return false
   return true

proc getH*(q: QRDecomposition): Matrix =
   ## Return the Householder vectors.
   ##
   ## ``return``: Lower trapezoidal matrix whose columns define the reflections.
   result = matrix(q.qr.m, q.qr.n) # sink here?!
   for i in 0 ..< q.qr.m:
      for j in 0 ..< q.qr.n:
         if i >= j:
            result[i, j] = q.qr[i, j]
         else:
            result[i, j] = 0.0

proc getR*(q: QRDecomposition): Matrix =
   ## Return the upper triangular factor.
   result = matrix(q.qr.n, q.qr.n)
   for i in 0 ..< q.qr.n:
      for j in 0 ..< q.qr.n:
         if i < j:
            result[i, j] = q.qr[i, j]
         elif i == j:
            result[i, j] = q.rDiag[i]
         else:
            result[i, j] = 0.0

proc getQ*(q: QRDecomposition): Matrix =
   ## Generate and return the (economy-sized) orthogonal factor.
   result = matrix(q.qr.m, q.qr.n)
   for k in countdown(q.qr.n - 1, 0):
      for i in 0 ..< q.qr.m:
         result[i, k] = 0.0
      result[k, k] = 1.0
      for j in k ..< q.qr.n:
         if q.qr[k, k] != 0.0:
            var s = 0.0
            for i in k ..< q.qr.m:
               s += q.qr[i, k] * result[i, j]
            s = -s / q.qr[k, k]
            for i in k ..< q.qr.m:
               result[i, j] += s * q.qr[i, k]

proc solve*(q: QRDecomposition, b: Matrix): Matrix =
   ## Least squares solution of ``A*X = B``
   ##
   ## - parameter ``b``: A Matrix with as many rows as A and any number of columns.
   ## - ``return``: ``X`` that minimizes the two norm of ``Q*R*X-B``
   assert(b.m == q.qr.m, "Matrix row dimensions must agree.")
   assert(q.isFullRank(), "Matrix is rank deficient.")
   # Copy right hand side
   var x = b
   # Compute Y = transpose(Q)*B
   for k in 0 ..< q.qr.n:
      for j in 0 ..< b.n:
         var s = 0.0
         for i in k ..< q.qr.m:
            s += q.qr[i, k] * x[i, j]
         s = -s / q.qr[k, k]
         for i in k ..< q.qr.m:
            x[i, j] += s * q.qr[i, k]
   # Solve R*X = Y
   for k in countdown(q.qr.n - 1, 0):
      for j in 0 ..< b.n:
         x[k, j] /= q.rDiag[k]
      for i in 0 ..< k:
         for j in 0 ..< b.n:
            x[i, j] -= x[k, j] * q.qr[i, k]
   x[0 ..< q.qr.n, 0 ..< b.n]
