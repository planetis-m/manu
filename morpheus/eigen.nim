## Eigenvalues and eigenvectors of a real matrix. 
##
## If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
## diagonal and the eigenvector matrix V is orthogonal.
## I.e. A = V.times(D.times(V.transpose())) and 
## V.times(V.transpose()) equals the identity matrix.
##
## If A is not symmetric, then the eigenvalue matrix D is block diagonal
## with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
## lambda + i * mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
## columns of V represent the eigenvectors in the sense that A*V = V*D,
## i.e. A.times(V) equals V.times(D).  The matrix V may be badly
## conditioned, or even singular, so the validity of the equation
## A = V*D*inverse(V) depends upon V.cond().
import "./matrix", math

type EigenvalueDecomposition* = object
   # Array for internal storage of eigenvectors.
   v: Matrix
   # Array for internal storage of nonsymmetric Hessenberg form.
   h: Matrix
   # Arrays for internal storage of eigenvalues.
   d, e: seq[float]
   # Working storage for nonsymmetric algorithm.
   ort: seq[float]
   # Row and column dimension (square matrix).
   n: int
   # Symmetry flag.
   issymmetric: bool

{.this: ei.}

proc tred2(ei: var EigenvalueDecomposition) =
   # Symmetric Householder reduction to tridiagonal form.
   #
   # This is derived from the Algol procedures tred2 by
   # Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   # Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   # Fortran subroutine in EISPACK.
   for j in 0 ..< n:
      d[j] = v[n - 1, j]
   # Householder reduction to tridiagonal form.
   for i in countdown(n - 1, 1):
      # Scale to avoid under/overflow.
      var scale = 0.0
      var h = 0.0
      for k in 0 ..< i:
         scale = scale + abs(d[k])
      if scale == 0.0:
         e[i] = d[i - 1]
         for j in 0 ..< i:
            d[j] = v[i - 1, j]
            v[i, j] = 0.0
            v[j, i] = 0.0
      else:
         # Generate Householder vector.
         for k in 0 ..< i:
            d[k] /= scale
            h += d[k] * d[k]
         var f = d[i - 1]
         var g = sqrt(h)
         if f > 0:
            g = -g
         e[i] = scale * g
         h = h - f * g
         d[i - 1] = f - g
         for j in 0 ..< i:
            e[j] = 0.0
         # Apply similarity transformation to remaining columns.
         for j in 0 ..< i:
            f = d[j]
            v[j, i] = f
            g = e[j] + v[j, j] * f
            for k in j + 1 .. i - 1:
               g += v[k, j] * d[k]
               e[k] += v[k, j] * f
            e[j] = g
         f = 0.0
         for j in 0 ..< i:
            e[j] /= h
            f += e[j] * d[j]
         let hh = f / (h + h)
         for j in 0 ..< i:
            e[j] -= hh * d[j]
         for j in 0 ..< i:
            f = d[j]
            g = e[j]
            for k in j .. i - 1:
               v[k, j] -= (f * e[k] + g * d[k])
            d[j] = v[i - 1, j]
            v[i, j] = 0.0
      d[i] = h
   # Accumulate transformations.
   for i in 0 .. n - 2:
      v[n - 1, i] = v[i, i]
      v[i, i] = 1.0
      var h = d[i + 1]
      if h != 0.0:
         for k in 0 .. i:
            d[k] = v[k, i + 1] / h
         for j in 0 .. i:
            var g = 0.0
            for k in 0 .. i:
               g += v[k, i + 1] * v[k, j]
            for k in 0 .. i:
               v[k, j] -= g * d[k]
      for k in 0 .. i:
         v[k, i + 1] = 0.0
   for j in 0 ..< n:
      d[j] = v[n - 1, j]
      v[n - 1, j] = 0.0
   v[n - 1, n - 1] = 1.0
   e[0] = 0.0

proc tql2(ei: var EigenvalueDecomposition) =
   # Symmetric tridiagonal QL algorithm.
   #
   # This is derived from the Algol procedures tql2, by
   # Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   # Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   # Fortran subroutine in EISPACK.
   for i in 1 ..< n:
      e[i - 1] = e[i]
   e[n - 1] = 0.0
   var f = 0.0
   var tst1 = 0.0
   let eps = pow(2.0, -52.0)
   for l in 0 ..< n:
      # Find small subdiagonal element
      tst1 = max(tst1, abs(d[l]) + abs(e[l]))
      var m = l
      while m < n:
         if abs(e[m]) <= eps * tst1:
            break
         m.inc
      # If m == l, d[l] is an eigenvalue,
      # otherwise, iterate.
      if m > l:
         var iter = 0
         while true:
            iter = iter + 1  # (Could check iteration count here.)
            # Compute implicit shift
            var g = d[l]
            var p = (d[l + 1] - g) / (2.0 * e[l])
            var r = hypot(p, 1.0)
            if p < 0:
               r = -r
            d[l] = e[l] / (p + r)
            d[l + 1] = e[l] * (p + r)
            let dl1 = d[l + 1]
            var h = g - d[l]
            for i in l+2 ..< n:
               d[i] -= h
            f = f + h
            # Implicit QL transformation.
            p = d[m]
            var c = 1.0
            var c2 = c
            var c3 = c
            let el1 = e[l + 1]
            var s = 0.0
            var s2 = 0.0
            for i in countdown(m - 1, l):
               c3 = c2
               c2 = c
               s2 = s
               g = c * e[i]
               h = c * p
               r = hypot(p, e[i])
               e[i + 1] = s * r
               s = e[i] / r
               c = p / r
               p = c * d[i] - s * g
               d[i + 1] = h + s * (c * g + s * d[i])
               # Accumulate transformation.
               for k in 0 ..< n:
                  h = v[k, i + 1]
                  v[k, i + 1] = s * v[k, i] + c * h
                  v[k, i] = c * v[k, i] - s * h
            p = -s * s2 * c3 * el1 * e[l] / dl1
            e[l] = s * p
            d[l] = c * p
            # Check for convergence.
            if abs(e[l]) <= eps * tst1: break
      d[l] = d[l] + f
      e[l] = 0.0
   # Sort eigenvalues and corresponding vectors.
   for i in 0 .. n - 2:
      var k = i
      var p = d[i]
      for j in i + 1 ..< n:
         if d[j] < p:
            k = j
            p = d[j]
      if k != i:
         d[k] = d[i]
         d[i] = p
         for j in 0 ..< n:
            p = v[j, i]
            v[j, i] = v[j, k]
            v[j, k] = p

proc orthes(ei: var EigenvalueDecomposition) =
   # Nonsymmetric reduction to Hessenberg form.
   #
   # This is derived from the Algol procedures orthes and ortran,
   # by Martin and Wilkinson, Handbook for Auto. Comp.,
   # Vol.ii-Linear Algebra, and the corresponding
   # Fortran subroutines in EISPACK.
   let low = 0
   let high = n - 1
   for m in low + 1 .. high - 1:
      # Scale column.
      var scale = 0.0
      for i in m .. high:
         scale = scale + abs(h[i, m - 1])
      if scale != 0.0:
         # Compute Householder transformation.
         var d = 0.0
         for i in countdown(high, m):
            ort[i] = h[i, m - 1] / scale
            d += ort[i] * ort[i]
         var g = sqrt(d)
         if ort[m] > 0:
            g = -g
         d = d - ort[m] * g
         ort[m] = ort[m] - g
         # Apply Householder similarity transformation
         # H = (I-u*u'/h)*H*(I-u*u')/h)
         for j in m ..< n:
            var f = 0.0
            for i in countdown(high, m):
               f += ort[i] * h[i, j]
            f = f / d
            for i in m .. high:
               h[i, j] -= f * ort[i]
         for i in 0 .. high:
            var f = 0.0
            for j in countdown(high, m):
               f += ort[j] * h[i, j]
            f = f / d
            for j in m .. high:
               h[i, j] -= f * ort[j]
         ort[m] = scale * ort[m]
         h[m, m - 1] = scale * g
   # Accumulate transformations (Algol's ortran).
   for i in 0 ..< n:
      for j in 0 ..< n:
         if i == j:
            v[i, j] = 1.0
         else:
            v[i, j] = 0.0
   for m in countdown(high - 1, low + 1):
      if h[m, m - 1] != 0.0:
         for i in m + 1 .. high:
            ort[i] = h[i, m - 1]
         for j in m .. high:
            var g = 0.0
            for i in m .. high:
               g += ort[i] * v[i, j]
            # Double division avoids possible underflow
            g = (g / ort[m]) / h[m, m - 1]
            for i in m .. high:
               v[i, j] += g * ort[i]

proc cdiv(xr, xi, yr, yi: float): tuple[re, im: float] =
   # Complex scalar division.
   var r, d: float
   if abs(yr) > abs(yi):
      r = yi / yr
      d = yr + r * yi
      result.re = (xr + r * xi) / d
      result.im = (xi - r * xr) / d
   else:
      r = yr / yi
      d = yi + r * yr
      result.re = (r * xr + xi) / d
      result.im = (r * xi - xr) / d

proc hqr2(ei: var EigenvalueDecomposition) =
   # Nonsymmetric reduction from Hessenberg to real Schur form.
   #
   # This is derived from the Algol procedure hqr2,
   # by Martin and Wilkinson, Handbook for Auto. Comp.,
   # Vol.ii-Linear Algebra, and the corresponding
   # Fortran subroutine in EISPACK.

   # Initialize   
   let nn = ei.n
   var n = nn - 1
   let low = 0
   let high = nn - 1
   let eps = pow(2.0, -52.0)
   var exshift = 0.0
   var p, q, r, s, z, t, w, x, y: float
   # Store roots isolated by balanc and compute matrix norm
   var norm = 0.0
   for i in 0 ..< nn:
      if i < low or i > high:
         d[i] = h[i, i]
         e[i] = 0.0
      for j in max(i - 1, 0) ..< nn:
         norm = norm + abs(h[i, j])
   # Outer loop over eigenvalue index
   var iter = 0
   while n >= low:
      # Look for single small sub-diagonal element
      var l = n
      while l > low:
         s = abs(h[l - 1, l - 1]) + abs(h[l, l])
         if s == 0.0:
            s = norm
         if abs(h[l, l - 1]) < eps * s:
            break
         l.dec
      # Check for convergence
      # One root found
      if l == n:
         h[n, n] = h[n, n] + exshift
         d[n] = h[n, n]
         e[n] = 0.0
         n.dec
         iter = 0
      # Two roots found
      elif l == n - 1:
         w = h[n, n - 1] * h[n - 1, n]
         p = (h[n - 1, n - 1] - h[n, n]) / 2.0
         q = p * p + w
         z = sqrt(abs(q))
         h[n, n] = h[n, n] + exshift
         h[n - 1, n - 1] = h[n - 1, n - 1] + exshift
         x = h[n, n]
         # Real pair
         if q >= 0:
            if p >= 0:
               z = p + z
            else:
               z = p - z
            d[n - 1] = x + z
            d[n] = d[n - 1]
            if z != 0.0:
               d[n] = x - w / z
            e[n - 1] = 0.0
            e[n] = 0.0
            x = h[n, n - 1]
            s = abs(x) + abs(z)
            p = x / s
            q = z / s
            r = sqrt(p * p + q * q)
            p = p / r
            q = q / r
            # Row modification
            for j in n - 1 ..< nn:
               z = h[n - 1, j]
               h[n - 1, j] = q * z + p * h[n, j]
               h[n, j] = q * h[n, j] - p * z
            # Column modification
            for i in 0 .. n:
               z = h[i, n - 1]
               h[i, n - 1] = q * z + p * h[i, n]
               h[i, n] = q * h[i, n] - p * z
            # Accumulate transformations
            for i in low .. high:
               z = v[i, n - 1]
               v[i, n - 1] = q * z + p * v[i, n]
               v[i, n] = q * v[i, n] - p * z
         # Complex pair
         else:
            d[n - 1] = x + p
            d[n] = x + p
            e[n - 1] = z
            e[n] = -z
         n = n - 2
         iter = 0
      # No convergence yet
      else:
         # Form shift
         x = h[n, n]
         y = 0.0
         w = 0.0
         if l < n:
            y = h[n - 1, n - 1]
            w = h[n, n - 1] * h[n - 1, n]
         # Wilkinson's original ad hoc shift
         if iter == 10:
            exshift += x
            for i in low .. n:
               h[i, i] -= x
            s = abs(h[n, n - 1]) + abs(h[n - 1, n - 2])
            y = 0.75 * s
            x = y
            w = -0.4375 * s * s
         # MATLAB's new ad hoc shift
         if iter == 30:
               s = (y - x) / 2.0
               s = s * s + w
               if s > 0:
                  s = sqrt(s)
                  if y < x:
                     s = -s
                  s = x - w / ((y - x) / 2.0 + s)
                  for i in low .. n:
                     h[i, i] -= s
                  exshift += s
                  w = 0.964
                  y = w
                  x = y
         iter.inc # (Could check iteration count here.)
         # Look for two consecutive small sub-diagonal elements
         var m = n - 2
         while m >= l:
            z = h[m, m]
            r = x - z
            s = y - z
            p = (r * s - w) / h[m + 1, m] + h[m, m + 1]
            q = h[m + 1, m + 1] - z - r - s
            r = h[m+2, m + 1]
            s = abs(p) + abs(q) + abs(r)
            p = p / s
            q = q / s
            r = r / s
            if m == l:
               break
            if abs(h[m, m - 1]) * (abs(q) + abs(r)) <
               eps * (abs(p) * (abs(h[m - 1, m - 1]) + abs(z) +
               abs(h[m + 1, m + 1]))):
                  break
            m.dec
         for i in m + 2 .. n:
            h[i, i - 2] = 0.0
            if i > m + 2:
               h[i, i - 3] = 0.0
         # Double QR step involving rows l:n and columns m:n
         for k in m .. n - 1:
            let notlast = k != n - 1
            if k != m:
               p = h[k, k - 1]
               q = h[k + 1, k - 1]
               r = if notlast: h[k+2, k - 1] else: 0.0
               x = abs(p) + abs(q) + abs(r)
               if x == 0.0:
                     continue
               p = p / x
               q = q / x
               r = r / x
            s = sqrt(p * p + q * q + r * r)
            if p < 0:
               s = -s
            if s != 0:
               if k != m:
                  h[k, k - 1] = -s * x
               elif l != m:
                  h[k, k - 1] = -h[k, k - 1]
               p = p + s
               x = p / s
               y = q / s
               z = r / s
               q = q / p
               r = r / p
               # Row modification
               for j in k ..< nn:
                  p = h[k, j] + q * h[k + 1, j]
                  if notlast:
                     p = p + r * h[k+2, j]
                     h[k+2, j] = h[k+2, j] - p * z
                  h[k, j] = h[k, j] - p * x
                  h[k + 1, j] = h[k + 1, j] - p * y
               # Column modification
               for i in 0 .. min(n, k + 3):
                  p = x * h[i, k] + y * h[i, k + 1]
                  if notlast:
                     p = p + z * h[i, k+2]
                     h[i, k+2] = h[i, k+2] - p * r

                  h[i, k] = h[i, k] - p
                  h[i, k + 1] = h[i, k + 1] - p * q
               # Accumulate transformations
               for i in low .. high:
                  p = x * v[i, k] + y * v[i, k + 1]
                  if notlast:
                     p = p + z * v[i, k+2]
                     v[i, k+2] = v[i, k+2] - p * r
                  v[i, k] = v[i, k] - p
                  v[i, k + 1] = v[i, k + 1] - p * q
            # (s != 0)
         # k loop
      # check convergence
   # while n >= low
   # Backsubstitute to find vectors of upper triangular form
   if norm == 0.0:
      return
   for n in countdown(nn - 1, 0):
      p = d[n]
      q = e[n]
      # Real vector
      if q == 0:
         var l = n
         h[n, n] = 1.0
         for i in countdown(n - 1, 0):
            w = h[i, i] - p
            r = 0.0
            for j in l .. n:
               r = r + h[i, j] * h[j, n]
            if e[i] < 0.0:
               z = w
               s = r
            else:
               l = i
               if e[i] == 0.0:
                  if w != 0.0:
                     h[i, n] = -r / w
                  else:
                     h[i, n] = -r / (eps * norm)
               # Solve real equations
               else:
                  x = h[i, i + 1]
                  y = h[i + 1, i]
                  q = (d[i] - p) * (d[i] - p) + e[i] * e[i]
                  t = (x * s - z * r) / q
                  h[i, n] = t
                  if abs(x) > abs(z):
                     h[i + 1, n] = (-r - w * t) / x
                  else:
                     h[i + 1, n] = (-s - y * t) / z
               # Overflow control
               t = abs(h[i, n])
               if (eps * t) * t > 1:
                  for j in i .. n:
                     h[j, n] = h[j, n] / t
      # Complex vector
      elif q < 0:
         var l = n - 1
         # Last vector component imaginary so matrix is triangular
         if abs(h[n, n - 1]) > abs(h[n - 1, n]):
            h[n - 1, n - 1] = q / h[n, n - 1]
            h[n - 1, n] = -(h[n, n] - p) / h[n, n - 1]
         else:
            let (cdivr, cdivi) = cdiv(0.0, -h[n - 1, n], h[n - 1, n - 1] - p, q)
            h[n - 1, n - 1] = cdivr
            h[n - 1, n] = cdivi
         h[n, n - 1] = 0.0
         h[n, n] = 1.0
         for i in countdown(n - 2, 0):
            var ra, sa, vr, vi: float
            ra = 0.0
            sa = 0.0
            for j in l .. n:
               ra = ra + h[i, j] * h[j, n - 1]
               sa = sa + h[i, j] * h[j, n]
            w = h[i, i] - p
            if e[i] < 0.0:
               z = w
               r = ra
               s = sa
            else:
               l = i
               if e[i] == 0:
                  let (cdivr, cdivi) = cdiv(-ra, -sa, w, q)
                  h[i, n - 1] = cdivr
                  h[i, n] = cdivi
               else:
                  # Solve complex equations
                  x = h[i, i + 1]
                  y = h[i + 1, i]
                  vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q
                  vi = (d[i] - p) * 2.0 * q
                  if vr == 0.0 and vi == 0.0:
                     vr = eps * norm * (abs(w) + abs(q) +
                           abs(x) + abs(y) + abs(z))
                  let (cdivr, cdivi) = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi)
                  h[i, n - 1] = cdivr
                  h[i, n] = cdivi
                  if abs(x) > (abs(z) + abs(q)):
                     h[i + 1, n - 1] = (-ra - w * h[i, n - 1] + q * h[i, n]) / x
                     h[i + 1, n] = (-sa - w * h[i, n] - q * h[i, n - 1]) / x
                  else:
                     let (cdivr, cdivi) = cdiv(-r - y * h[i, n - 1], -s - y * h[i, n], z, q)
                     h[i + 1, n - 1] = cdivr
                     h[i + 1, n] = cdivi
               # Overflow control
               t = max(abs(h[i, n - 1]), abs(h[i, n]))
               if (eps * t) * t > 1:
                  for j in i .. n:
                     h[j, n - 1] = h[j, n - 1] / t
                     h[j, n] = h[j, n] / t
   # Vectors of isolated roots
   for i in 0 ..< nn:
      if i < low or i > high:
         for j in i ..< nn:
            v[i, j] = h[i, j]
   # Back transformation to get eigenvectors of original matrix
   for j in countdown(nn - 1, low):
      for i in low .. high:
         z = 0.0
         for k in low .. min(j, high):
            z = z + v[i, k] * h[k, j]
         v[i, j] = z

proc eig*(a: Matrix): EigenvalueDecomposition =
   ## Check for symmetry, then construct the eigenvalue decomposition
   ##
   ## - ``return``: Structure to access D and V.
   ## - parameter ``a``: Square matrix
   let n = a.n
   result.n = n
   result.d = newSeq[float](n)
   result.e = newSeq[float](n)
   var issymmetric = true
   for j in 0 ..< n:
      if not issymmetric:
         break
      for i in 0 ..< n:
         if not issymmetric:
            break
         issymmetric = a[i, j] == a[j, i]
   result.issymmetric = issymmetric
   if issymmetric:
      result.v = a
      # Tridiagonalize.
      result.tred2()
      # Diagonalize.
      result.tql2()
   else:
      result.h = a
      result.v = matrix(n, n)
      result.ort = newSeq[float](n)
      # Reduce to Hessenberg form.
      result.orthes()
      # Reduce Hessenberg to real Schur form.
      result.hqr2()

proc getV*(ei: EigenvalueDecomposition): Matrix {.inline.} =
   ## Return the eigenvector matrix
   v

proc getRealEigenvalues*(ei: EigenvalueDecomposition): seq[float] {.inline.} =
   ## Return the real parts of the eigenvalues
   ##
   ## ``return``: real(diag(D))
   d

proc getImagEigenvalues*(ei: EigenvalueDecomposition): seq[float] {.inline.} =
   ## Return the imaginary parts of the eigenvalues
   ##
   ## ``return``: imag(diag(D))
   e

proc getD*(ei: EigenvalueDecomposition): Matrix =
   ## Return the block diagonal eigenvalue matrix
   result = matrix(n, n)
   for i in 0 ..< n:
      # for j in 0 ..< n:
      #    result[i, j] = 0.0
      result[i, i] = d[i]
      if e[i] > 0:
         result[i, i + 1] = e[i]
      elif e[i] < 0:
         result[i, i - 1] = e[i]
