## Singular Value Decomposition.
##
## For an m-by-n matrix A with m >= n, the singular value decomposition is
## an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
## an n-by-n orthogonal matrix V so that A = U*S*V'.
##
## The singular values, sigma[k] = S[k][k], are ordered so that
## sigma[0] >= sigma[1] >= ... >= sigma[n-1].
##
## The singular value decompostion always exists, so the constructor will
## never fail.  The matrix condition number and the effective numerical
## rank can be computed from this decomposition.
import matrix

type SingularValueDecomposition* = object
   # Arrays for internal storage of U and V.
   u, v: seq[seq[float]]
   # Array for internal storage of singular values.
   s: seq[float]
   # Row and column dimensions.
   m, n: int

proc svd*(m: Matrix): SingularValueDecomposition =
   ## Construct the singular value decomposition
   ##  Structure to access U, S and V.

   # Derived from LINPACK code.
   # Initialize.
   var a = m.data
   result.m = m.m
   result.n = m.n

   let nu = min(m, n)
   result.s = newSeq[float](min(m + 1, n))
   result.u = newSeq[float](m)[nu]
   result.v = newSeq[float](n)[n]
   var e = newSeq[float](n)
   var work = newSeq[float](m)
   var wantu = true
   var wantv = true

   # Reduce A to bidiagonal form, storing the diagonal elements
   # in s and the super-diagonal elements in e.
   let nct = min(m - 1, n)
   let nrt = max(0, min(n - 2, m))
   for k in 0 ..< max(nct, nrt):
      if k < nct:
         # Compute the transformation for the k-th column and
         # place the k-th diagonal in s[k].
         # Compute 2-norm of k-th column without under/overflow.
         s[k] = 0
         for i in k ..< m:
            s[k] = hypot(s[k], a[i][k])
         if s[k] != 0.0:
            if a[k][k] < 0.0:
               s[k] = -s[k]
            for i in k ..< m:
               a[i][k] /= s[k]
            a[k][k] += 1.0
         s[k] = -s[k]
      for j in k + 1 ..< n:
         if k < nct and s[k] != 0.0:
            # Apply the transformation.
            var t = 0.0
            for i in k ..< m:
               t += a[i][k] * a[i][j]
            t = -t / a[k][k]
            for i in k ..< m:
               a[i][j] += t * a[i][k]
         # Place the k-th row of A into e for the
         # subsequent calculation of the row transformation.
         e[j] = a[k][j]
      if wantu and k < nct:
         # Place the transformation in U for subsequent back
         # multiplication.
         for i in k ..< m:
            result.u[i][k] = a[i][k]
      if k < nrt:
         # Compute the k-th row transformation and place the
         # k-th super-diagonal in e[k].
         # Compute 2-norm without under/overflow.
         e[k] = 0
         for i in k + 1 ..< n:
            e[k] = hypot(e[k], e[i])
         if e[k] != 0.0:
            if e[k + 1] < 0.0:
               e[k] = -e[k]
            for i in k + 1 ..< n:
               e[i] /= e[k]
            e[k + 1] += 1.0
         e[k] = -e[k]
         if k + 1 < m and e[k] != 0.0:
            # Apply the transformation.
            for i in k + 1 ..< m:
               work[i] = 0.0
            for j in k + 1 ..< n:
               for i in k + 1 ..< m:
                  work[i] += e[j] * a[i][j]
            for j in k + 1 ..< n:
               var t = -e[j] / e[k + 1]
               for i in k + 1 ..< m:
                  a[i][j] += t * work[i]
         if wantv:
            # Place the transformation in V for subsequent
            # back multiplication.
            for i in k + 1 ..< n:
               result.v[i][k] = e[i]

   # Set up the final bidiagonal matrix or order p.
   let p = min(n, m + 1)
   if nct < n:
      s[nct] = a[nct][nct]
   if m < p:
      s[p - 1] = 0.0
   if nrt + 1 < p:
      e[nrt] = a[nrt][p - 1]
   e[p - 1] = 0.0

   # If required, generate U.
   if wantu:
      for j in nct ..< nu:
         for i in 0 ..< m:
            result.u[i][j] = 0.0
         result.u[j][j] = 1.0
      for k in countdown(nct - 1, 0):
         if s[k] != 0.0:
            for j in k + 1 ..< nu:
               var t = 0.0
               for i in k ..< m:
                  t += result.u[i][k] * result.u[i][j]
               t = -t / result.u[k][k]
               for i in k ..< m:
                  result.u[i][j] += t * result.u[i][k]
            for i in k ..< m:
               result.u[i][k] = -result.u[i][k]
            result.u[k][k] = 1.0 + result.u[k][k]
            for i in 0 .. k - 2:
               result.u[i][k] = 0.0
         else:
            for i in 0 ..< m:
               result.u[i][k] = 0.0
            result.u[k][k] = 1.0

   # If required, generate V.
   if wantv:
      for k in countdown(n - 1, 0):
         if k < nrt and e[k] != 0.0:
            for j in k + 1 ..< nu:
               var t = 0.0
               for i in k + 1 ..< n:
                  t += result.v[i][k] * result.v[i][j]
               t = -t / result.v[k + 1][k]
               for i in k + 1 ..< n:
                  result.v[i][j] += t * result.v[i][k]
         for i in 0 ..< n:
            result.v[i][k] = 0.0
         result.v[k][k] = 1.0

   # Main iteration loop for the singular values.
   let pp = p - 1
   var iter = 0
   let eps = pow(2.0, -52.0)
   let tiny = pow(2.0, -966.0)
   while p > 0:
      var k, kase: int
      # Here is where a test for too many iterations would go.

      # This section of the program inspects for
      # negligible elements in the s and e arrays.  On
      # completion the variables kase and k are set as follows.

      # kase = 1     if s(p) and e[k-1] are negligible and k<p
      # kase = 2     if s(k) is negligible and k<p
      # kase = 3     if e[k-1] is negligible, k<p, and
      #              s(k), ..., s(p) are not negligible (qr step).
      # kase = 4     if e(p-1) is negligible (convergence).
      for k in countdown(p - 2, -1):
         if k == -1
            break
         if abs(e[k]) <=
               tiny + eps * (abs(s[k]) + abs(s[k + 1])):
            e[k] = 0.0
            break
      if k == p - 2:
         kase = 4
      else:
         var ks: int
         for ks in countdown(p - 1, k):
            if ks == k:
               break
            let t = (ks != p ? abs(e[ks]) : 0.0) + 
                        (ks != k+1 ? abs(e[ks-1]) : 0.0)
            if abs(s[ks]) <= tiny + eps * t:
               result.s[ks] = 0.0
               break
         if ks == k:
            kase = 3
         elif ks == p - 1:
            kase = 1
         else:
            kase = 2
            k = ks
      k.inc

      # Perform the task indicated by kase.

      case kase
         # Deflate negligible s(p).
         of 1:
            double f = e[p-2];
            e[p-2] = 0.0;
            for (int j = p-2; j >= k; j--) {
               double t = Maths.hypot(s[j],f);
               double cs = s[j]/t;
               double sn = f/t;
               s[j] = t;
               if (j != k) {
                  f = -sn*e[j-1];
                  e[j-1] = cs*e[j-1];
               }
               if (wantv) {
                  for (int i = 0; i < n; i++) {
                     t = cs*V[i][j] + sn*V[i][p-1];
                     V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                     V[i][j] = t
         # Split at negligible s(k).
         of 2:
            double f = e[k-1];
            e[k-1] = 0.0;
            for (int j = k; j < p; j++) {
               double t = Maths.hypot(s[j],f);
               double cs = s[j]/t;
               double sn = f/t;
               s[j] = t;
               f = -sn*e[j];
               e[j] = cs*e[j];
               if (wantu) {
                  for (int i = 0; i < m; i++) {
                     t = cs*U[i][j] + sn*U[i][k-1];
                     U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                     U[i][j] = t
         # Perform one qr step.
         of 3:
            # Calculate the shift.
            double scale = Math.max(Math.max(Math.max(Math.max(
                     Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])), 
                     Math.abs(s[k])),Math.abs(e[k]));
            double sp = s[p-1]/scale;
            double spm1 = s[p-2]/scale;
            double epm1 = e[p-2]/scale;
            double sk = s[k]/scale;
            double ek = e[k]/scale;
            double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
            double c = (sp*epm1)*(sp*epm1);
            double shift = 0.0;
            if b != 0.0 | c != 0.0:
               shift = sqrt(b * b + c)
               if b < 0.0:
                  shift = -shift
               shift = c / (b + shift)
            double f = (sk + sp) * (sk - sp) + shift
            double g = sk * ek

            # Chase zeros.
            for (int j = k; j < p-1; j++) {
               double t = Maths.hypot(f,g);
               double cs = f/t;
               double sn = g/t;
               if (j != k) {
                  e[j-1] = t;
               }
               f = cs*s[j] + sn*e[j];
               e[j] = cs*e[j] - sn*s[j];
               g = sn*s[j+1];
               s[j+1] = cs*s[j+1];
               if (wantv) {
                  for (int i = 0; i < n; i++) {
                     t = cs*V[i][j] + sn*V[i][j+1];
                     V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                     V[i][j] = t;
                  }
               }
               t = Maths.hypot(f,g);
               cs = f/t;
               sn = g/t;
               s[j] = t;
               f = cs*e[j] + sn*s[j+1];
               s[j+1] = -sn*e[j] + cs*s[j+1];
               g = sn*e[j+1];
               e[j+1] = cs*e[j+1];
               if (wantu && (j < m-1)) {
                  for (int i = 0; i < m; i++) {
                     t = cs*U[i][j] + sn*U[i][j+1];
                     U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                     U[i][j] = t;
                  }
               }
            }
            e[p-2] = f
            iter = iter + 1
         # Convergence.
         of 4:
            # Make the singular values positive.
            if s[k] <= 0.0:
               s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
               if (wantv) {
                  for (int i = 0; i <= pp; i++) {
                     V[i][k] = -V[i][k]
            # Order the singular values.
            while k < pp:
               if s[k] >= s[k + 1]:
                  break
               double t = s[k];
               s[k] = s[k+1];
               s[k+1] = t;
               if (wantv && (k < n-1)) {
                  for (int i = 0; i < n; i++) {
                     t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t
               if (wantu && (k < m-1)) {
                  for (int i = 0; i < m; i++) {
                     t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t
               k++
            iter = 0
            p--

proc getU*(s: SingularValueDecomposition): Matrix =
   ## Return the left singular vectors
   return new Matrix(U,m,Math.min(m+1,n))

proc getV*(s: SingularValueDecomposition): Matrix =
   ## Return the right singular vectors
   return new Matrix(V,n,n)

proc getSingularValues*(s: SingularValueDecomposition): seq[float] =
   ## Return the one-dimensional array of singular values
   ## returns diagonal of S.
   result = s.s

proc getS*(s: SingularValueDecomposition): Matrix =
   ## Return the diagonal matrix of singular values
   result.m = n
   result.n = n
   newData()
   for i in 0 ..< s.n:
      for j in 0 ..< s.n:
         result.data[i][j] = 0.0
      result.data[i][i] = s.s[i]

proc norm2*(s: SingularValueDecomposition): float =
   ## Two norm
   ## returns max(S)
   result = s[0]

proc cond*(s: SingularValueDecomposition): float =
   ## Two norm condition number
   ## return max(S)/min(S)
   s.s[0] / s.s[min(m, n) - 1]
   
proc rank*(s: SingularValueDecomposition): int =
   ## Effective numerical matrix rank
   ## returns Number of nonnegligible singular values.
   let eps = pow(2.0, -52.0)
   let tol = max(m, n) * s.s[0] * eps
   for i in 0 ..< s.len:
      if s.s[i] > tol:
         result.inc
