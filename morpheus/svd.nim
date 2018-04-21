## Singular Value Decomposition.
##
## For an m-by-n matrix A with m >= n, the singular value decomposition is
## an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
## an n-by-n orthogonal matrix V so that A = U*S*V'.
##
## The singular values, sigma[k] = S[k, k], are ordered so that
## sigma[0] >= sigma[1] >= ... >= sigma[n-1].
##
## The singular value decompostion always exists, so the constructor will
## never fail.  The matrix condition number and the effective numerical
## rank can be computed from this decomposition.
import math
import "./matrix"

type SingularValueDecomposition* = object
   # Arrays for internal storage of U and V.
   u, v: Matrix
   # Array for internal storage of singular values.
   sv: seq[float]

{.this: s.}

proc svd*(a: Matrix): SingularValueDecomposition =
   ## Construct the singular value decomposition
   ## Structure to access U, S and V.
   # Derived from LINPACK code.
   # Initialize.
   var a = a
   let m = a.m
   let n = a.n
   let nu = min(m, n)
   result.sv = newSeq[float](min(m + 1, n))
   result.u = matrix(m, nu)
   result.v = matrix(n, n)
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
         # place the k-th diagonal in sv[k].
         # Compute 2-norm of k-th column without under/overflow.
         result.sv[k] = 0
         for i in k ..< m:
            result.sv[k] = hypot(result.sv[k], a[i, k])
         if result.sv[k] != 0.0:
            if a[k, k] < 0.0:
               result.sv[k] = -result.sv[k]
            for i in k ..< m:
               a[i, k] /= result.sv[k]
            a[k, k] += 1.0
         result.sv[k] = -result.sv[k]
      for j in k + 1 ..< n:
         if k < nct and result.sv[k] != 0.0:
            # Apply the transformation.
            var t = 0.0
            for i in k ..< m:
               t += a[i, k] * a[i, j]
            t = -t / a[k, k]
            for i in k ..< m:
               a[i, j] += t * a[i, k]
         # Place the k-th row of A into e for the
         # subsequent calculation of the row transformation.
         e[j] = a[k, j]
      if wantu and k < nct:
         # Place the transformation in U for subsequent back
         # multiplication.
         for i in k ..< m:
            result.u[i, k] = a[i, k]
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
                  work[i] += e[j] * a[i, j]
            for j in k + 1 ..< n:
               var t = -e[j] / e[k + 1]
               for i in k + 1 ..< m:
                  a[i, j] += t * work[i]
         if wantv:
            # Place the transformation in V for subsequent
            # back multiplication.
            for i in k + 1 ..< n:
               result.v[i, k] = e[i]
   # Set up the final bidiagonal matrix or order p.
   var p = min(n, m + 1)
   if nct < n:
      result.sv[nct] = a[nct, nct]
   if m < p:
      result.sv[p - 1] = 0.0
   if nrt + 1 < p:
      e[nrt] = a[nrt, p - 1]
   e[p - 1] = 0.0
   # If required, generate U.
   if wantu:
      for j in nct ..< nu:
         for i in 0 ..< m:
            result.u[i, j] = 0.0
         result.u[j, j] = 1.0
      for k in countdown(nct - 1, 0):
         if result.sv[k] != 0.0:
            for j in k + 1 ..< nu:
               var t = 0.0
               for i in k ..< m:
                  t += result.u[i, k] * result.u[i, j]
               t = -t / result.u[k, k]
               for i in k ..< m:
                  result.u[i, j] += t * result.u[i, k]
            for i in k ..< m:
               result.u[i, k] = -result.u[i, k]
            result.u[k, k] = 1.0 + result.u[k, k]
            for i in 0 .. k - 2:
               result.u[i, k] = 0.0
         else:
            for i in 0 ..< m:
               result.u[i, k] = 0.0
            result.u[k, k] = 1.0
   # If required, generate V.
   if wantv:
      for k in countdown(n - 1, 0):
         if k < nrt and e[k] != 0.0:
            for j in k + 1 ..< nu:
               var t = 0.0
               for i in k + 1 ..< n:
                  t += result.v[i, k] * result.v[i, j]
               t = -t / result.v[k + 1, k]
               for i in k + 1 ..< n:
                  result.v[i, j] += t * result.v[i, k]
         for i in 0 ..< n:
            result.v[i, k] = 0.0
         result.v[k, k] = 1.0
   # Main iteration loop for the singular values.
   let pp = p - 1
   var iter = 0
   let eps = pow(2.0, -52.0)
   let tiny = pow(2.0, -966.0)
   while p > 0:
      var kase = 0
      # Here is where a test for too many iterations would go.

      # This section of the program inspects for
      # negligible elements in the s and e arrays.  On
      # completion the variables kase and k are set as follows.

      # kase = 1     if s(p) and e[k-1] are negligible and k<p
      # kase = 2     if s(k) is negligible and k<p
      # kase = 3     if e[k-1] is negligible, k<p, and
      #              s(k), ..., s(p) are not negligible (qr step).
      # kase = 4     if e(p-1) is negligible (convergence).
      var k = p - 2
      while k >= 0:
         if abs(e[k]) <=
               tiny + eps * (abs(result.sv[k]) + abs(result.sv[k + 1])):
            e[k] = 0.0
            break
         k.dec
      if k == p - 2:
         kase = 4
      else:
         var ks = p - 1
         while ks > k:
            var t = 0.0
            if ks != p: t += abs(e[ks])
            if ks != k + 1: t += abs(e[ks-1])
            if abs(result.sv[ks]) <= tiny + eps * t:
               result.sv[ks] = 0.0
               break
            ks.dec
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
         var f = e[p - 2]
         e[p-2] = 0.0
         for j in countdown(p-2, k):
            var t = hypot(result.sv[j], f)
            let cs = result.sv[j] / t
            let sn = f / t
            result.sv[j] = t
            if j != k:
               f = -sn * e[j - 1]
               e[j - 1] = cs * e[j - 1]
            if wantv:
               for i in 0 ..< n:
                  t = cs * result.v[i, j] + sn * result.v[i, p - 1]
                  result.v[i, p - 1] = -sn * result.v[i, j] + cs * result.v[i, p - 1]
                  result.v[i, j] = t
      # Split at negligible s(k).
      of 2:
         var f = e[k - 1]
         e[k - 1] = 0.0
         for j in k ..< p:
            var t = hypot(result.sv[j], f)
            let cs = result.sv[j] / t
            let sn = f / t
            result.sv[j] = t
            f = -sn * e[j]
            e[j] = cs * e[j]
            if wantu:
               for i in 0 ..< m:
                  t = cs * result.u[i, j] + sn * result.u[i, k - 1]
                  result.u[i, k - 1] = -sn * result.u[i, j] + cs * result.u[i, k - 1]
                  result.u[i, j] = t
      # Perform one qr step.
      of 3:
         # Calculate the shift.
         let scale = max(max(max(max(
                  abs(result.sv[p - 1]), abs(result.sv[p - 2])), abs(e[p - 2])), 
                  abs(result.sv[k])), abs(e[k]))
         let sp = result.sv[p-1] / scale
         let spm1 = result.sv[p-2] / scale
         let epm1 = e[p-2] / scale
         let sk = result.sv[k] / scale
         let ek = e[k] / scale
         let b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0
         let c = (sp * epm1) * (sp * epm1)
         var shift = 0.0
         if b != 0.0 or c != 0.0:
            shift = sqrt(b * b + c)
            if b < 0.0:
               shift = -shift
            shift = c / (b + shift)
         var f = (sk + sp) * (sk - sp) + shift
         var g = sk * ek
         # Chase zeros.
         for j in k ..< p - 1:
            var t = hypot(f, g)
            var cs = f / t
            var sn = g / t
            if j != k:
               e[j-1] = t
            f = cs * result.sv[j] + sn * e[j]
            e[j] = cs * e[j] - sn * result.sv[j]
            g = sn * result.sv[j + 1]
            result.sv[j + 1] = cs * result.sv[j + 1]
            if wantv:
               for i in 0 ..< n:
                  t = cs * result.v[i, j] + sn * result.v[i, j + 1]
                  result.v[i, j + 1] = -sn * result.v[i, j] + cs * result.v[i, j + 1]
                  result.v[i, j] = t
            t = hypot(f, g)
            cs = f / t
            sn = g / t
            result.sv[j] = t
            f = cs * e[j] + sn * result.sv[j + 1]
            result.sv[j + 1] = -sn * e[j] + cs * result.sv[j + 1]
            g = sn * e[j + 1]
            e[j + 1] = cs * e[j + 1]
            if wantu and j < m - 1:
               for i in 0 ..< m:
                  t = cs * result.u[i, j] + sn * result.u[i, j + 1]
                  result.u[i, j + 1] = -sn * result.u[i, j] + cs * result.u[i, j + 1]
                  result.u[i, j] = t
         e[p - 2] = f
         iter.inc
      # Convergence.
      of 4:
         # Make the singular values positive.
         if result.sv[k] <= 0.0:
            result.sv[k] = if result.sv[k] < 0.0: -result.sv[k] else: 0.0
            if wantv:
               for i in 0 .. pp:
                  result.v[i, k] = -result.v[i, k]
         # Order the singular values.
         while k < pp:
            if result.sv[k] >= result.sv[k + 1]:
               break
            swap(result.sv[k+1], result.sv[k])
            if wantv and k < n - 1:
               for i in 0 ..< n:
                  swap(result.v[i, k + 1], result.v[i, k])
            if wantu and k < m - 1:
               for i in 0 ..< m:
                  swap(result.u[i, k + 1], result.u[i, k])
            k.inc
         iter = 0
         p.dec
      else: discard

proc getU*(s: SingularValueDecomposition): Matrix =
   ## Return the left singular vectors
   u

proc getV*(s: SingularValueDecomposition): Matrix =
   ## Return the right singular vectors
   v

proc getSingularValues*(s: SingularValueDecomposition): seq[float] =
   ## Return the one-dimensional array of singular values.
   ## return: diagonal of S
   sv

proc getS*(s: SingularValueDecomposition): Matrix =
   ## Return the diagonal matrix of singular values.
   result = matrix(v.m, v.m)
   for i in 0 ..< v.m:
      # for j in 0 ..< n:
      #    result.data[i, j] = 0.0
      result[i, i] = sv[i]

proc norm2*(s: SingularValueDecomposition): float =
   ## Two norm.
   ## return: max(S)
   sv[0]

proc cond*(s: SingularValueDecomposition): float =
   ## Two norm condition number.
   ## return: max(S)/min(S)
   sv[0] / sv[^1]

proc rank*(s: SingularValueDecomposition): int =
   ## Effective numerical matrix rank.
   ## return: Number of nonnegligible singular values.
   let eps = pow(2.0, -52.0)
   let tol = float(max(u.m, v.m)) * sv[0] * eps
   for d in sv:
      if d > tol:
         result.inc
