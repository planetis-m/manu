## Eigenvalues and eigenvectors of a real matrix
## =============================================
##
## If A is symmetric, then ``A = V*D*V'`` where the eigenvalue matrix D is
## diagonal and the eigenvector matrix V is orthogonal.
## I.e. ``A = V*D*V.transpose()`` and
## ``V*V.transpose()`` equals the identity matrix.
##
## If A is not symmetric, then the eigenvalue matrix D is block diagonal
## with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
## ``lambda + i * mu``, in 2-by-2 blocks, ``[lambda, mu; -mu, lambda]``.  The
## columns of V represent the eigenvectors in the sense that ``A*V = V*D``,
## i.e. ``A*V`` equals ``V*D``. The matrix V may be badly
## conditioned, or even singular, so the validity of the equation
## ``A = V*D*inverse(V)`` depends upon `V.cond()`.
import matrix, std/[math, fenv]

type
  EigenvalueDecomposition*[T] = object
    v: Matrix[T] # Array for internal storage of eigenvectors.
    h: Matrix[T] # Array for internal storage of nonsymmetric Hessenberg form.
    d, e: seq[T] # Arrays for internal storage of eigenvalues.
    ort: seq[T]  # Working storage for nonsymmetric algorithm.

proc tred2[T](ei: var EigenvalueDecomposition[T]) =
  # Symmetric Householder reduction to tridiagonal form.
  #
  # This is derived from the Algol procedures tred2 by
  # Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  # Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  # Fortran subroutine in EISPACK.
  let n = ei.v.n
  for j in 0 ..< n:
    ei.d[j] = ei.v[n - 1, j]
  # Householder reduction to tridiagonal form.
  for i in countdown(n - 1, 1):
    # Scale to avoid under/overflow.
    var scale = T(0)
    var h = T(0)
    for k in 0 ..< i:
      scale = scale + abs(ei.d[k])
    if scale == 0:
      ei.e[i] = ei.d[i - 1]
      for j in 0 ..< i:
        ei.d[j] = ei.v[i - 1, j]
        ei.v[i, j] = T(0)
        ei.v[j, i] = T(0)
    else:
      # Generate Householder vector.
      for k in 0 ..< i:
        ei.d[k] /= scale
        h += ei.d[k] * ei.d[k]
      var f = ei.d[i - 1]
      var g = sqrt(h)
      if f > 0:
        g = -g
      ei.e[i] = scale * g
      h = h - f * g
      ei.d[i - 1] = f - g
      for j in 0 ..< i:
        ei.e[j] = T(0)
      # Apply similarity transformation to remaining columns.
      for j in 0 ..< i:
        f = ei.d[j]
        ei.v[j, i] = f
        g = ei.e[j] + ei.v[j, j] * f
        for k in j + 1 .. i - 1:
          g += ei.v[k, j] * ei.d[k]
          ei.e[k] += ei.v[k, j] * f
        ei.e[j] = g
      f = T(0)
      for j in 0 ..< i:
        ei.e[j] /= h
        f += ei.e[j] * ei.d[j]
      let hh = f / (h + h)
      for j in 0 ..< i:
        ei.e[j] -= hh * ei.d[j]
      for j in 0 ..< i:
        f = ei.d[j]
        g = ei.e[j]
        for k in j .. i - 1:
          ei.v[k, j] -= (f * ei.e[k] + g * ei.d[k])
        ei.d[j] = ei.v[i - 1, j]
        ei.v[i, j] = T(0)
    ei.d[i] = h
  # Accumulate transformations.
  for i in 0 .. n - 2:
    ei.v[n - 1, i] = ei.v[i, i]
    ei.v[i, i] = T(1)
    var h = ei.d[i + 1]
    if h != 0:
      for k in 0 .. i:
        ei.d[k] = ei.v[k, i + 1] / h
      for j in 0 .. i:
        var g = T(0)
        for k in 0 .. i:
          g += ei.v[k, i + 1] * ei.v[k, j]
        for k in 0 .. i:
          ei.v[k, j] -= g * ei.d[k]
    for k in 0 .. i:
      ei.v[k, i + 1] = T(0)
  for j in 0 ..< n:
    ei.d[j] = ei.v[n - 1, j]
    ei.v[n - 1, j] = T(0)
  ei.v[n - 1, n - 1] = T(1)
  ei.e[0] = T(0)

proc tql2[T](ei: var EigenvalueDecomposition[T]) =
  # Symmetric tridiagonal QL algorithm.
  #
  # This is derived from the Algol procedures tql2, by
  # Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  # Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  # Fortran subroutine in EISPACK.
  let n = ei.v.n
  for i in 1 ..< n:
    ei.e[i - 1] = ei.e[i]
  ei.e[n - 1] = T(0)
  var f = T(0)
  var tst1 = T(0)
  let eps = epsilon(T)
  for l in 0 ..< n:
    # Find small subdiagonal element
    tst1 = max(tst1, abs(ei.d[l]) + abs(ei.e[l]))
    var m = l
    while m < n:
      if abs(ei.e[m]) <= eps * tst1:
        break
      m.inc
    # If m == l, d[l] is an eigenvalue,
    # otherwise, iterate.
    if m > l:
      var iter = 0
      while true:
        iter.inc # (Could check iteration count here.)
            # Compute implicit shift
        var g = ei.d[l]
        var p = (ei.d[l + 1] - g) / (T(2) * ei.e[l])
        var r = hypot(p, T(1))
        if p < 0:
          r = -r
        ei.d[l] = ei.e[l] / (p + r)
        ei.d[l + 1] = ei.e[l] * (p + r)
        let dl1 = ei.d[l + 1]
        var h = g - ei.d[l]
        for i in l+2 ..< n:
          ei.d[i] -= h
        f = f + h
        # Implicit QL transformation.
        p = ei.d[m]
        var c = T(1)
        var c2 = c
        var c3 = c
        let el1 = ei.e[l + 1]
        var s = T(0)
        var s2 = T(0)
        for i in countdown(m - 1, l):
          c3 = c2
          c2 = c
          s2 = s
          g = c * ei.e[i]
          h = c * p
          r = hypot(p, ei.e[i])
          ei.e[i + 1] = s * r
          s = ei.e[i] / r
          c = p / r
          p = c * ei.d[i] - s * g
          ei.d[i + 1] = h + s * (c * g + s * ei.d[i])
          # Accumulate transformation.
          for k in 0 ..< n:
            h = ei.v[k, i + 1]
            ei.v[k, i + 1] = s * ei.v[k, i] + c * h
            ei.v[k, i] = c * ei.v[k, i] - s * h
        p = -s * s2 * c3 * el1 * ei.e[l] / dl1
        ei.e[l] = s * p
        ei.d[l] = c * p
        # Check for convergence.
        if abs(ei.e[l]) <= eps * tst1: break
    ei.d[l] = ei.d[l] + f
    ei.e[l] = T(0)
  # Sort eigenvalues and corresponding vectors.
  for i in 0 .. n - 2:
    var k = i
    var p = ei.d[i]
    for j in i + 1 ..< n:
      if ei.d[j] < p:
        k = j
        p = ei.d[j]
    if k != i:
      ei.d[k] = ei.d[i]
      ei.d[i] = p
      for j in 0 ..< n:
        p = ei.v[j, i]
        ei.v[j, i] = ei.v[j, k]
        ei.v[j, k] = p

proc orthes[T](ei: var EigenvalueDecomposition[T]) =
  # Nonsymmetric reduction to Hessenberg form.
  #
  # This is derived from the Algol procedures orthes and ortran,
  # by Martin and Wilkinson, Handbook for Auto. Comp.,
  # Vol.ii-Linear Algebra, and the corresponding
  # Fortran subroutines in EISPACK.
  let n = ei.v.n
  let low = 0
  let high = n - 1
  for m in low + 1 .. high - 1:
    # Scale column.
    var scale = T(0)
    for i in m .. high:
      scale = scale + abs(ei.h[i, m - 1])
    if scale != 0:
      # Compute Householder transformation.
      var d = T(0)
      for i in countdown(high, m):
        ei.ort[i] = ei.h[i, m - 1] / scale
        d += ei.ort[i] * ei.ort[i]
      var g = sqrt(d)
      if ei.ort[m] > 0:
        g = -g
      d = d - ei.ort[m] * g
      ei.ort[m] = ei.ort[m] - g
      # Apply Householder similarity transformation
      # H = (I-u*u'/h)*H*(I-u*u')/h)
      for j in m ..< n:
        var f = T(0)
        for i in countdown(high, m):
          f += ei.ort[i] * ei.h[i, j]
        f = f / d
        for i in m .. high:
          ei.h[i, j] -= f * ei.ort[i]
      for i in 0 .. high:
        var f = T(0)
        for j in countdown(high, m):
          f += ei.ort[j] * ei.h[i, j]
        f = f / d
        for j in m .. high:
          ei.h[i, j] -= f * ei.ort[j]
      ei.ort[m] = scale * ei.ort[m]
      ei.h[m, m - 1] = scale * g
  # Accumulate transformations (Algol's ortran).
  for i in 0 ..< n:
    for j in 0 ..< n:
      if i == j:
        ei.v[i, j] = T(1)
      else:
        ei.v[i, j] = T(0)
  for m in countdown(high - 1, low + 1):
    if ei.h[m, m - 1] != 0:
      for i in m + 1 .. high:
        ei.ort[i] = ei.h[i, m - 1]
      for j in m .. high:
        var g = T(0)
        for i in m .. high:
          g += ei.ort[i] * ei.v[i, j]
        # Double division avoids possible underflow
        g = (g / ei.ort[m]) / ei.h[m, m - 1]
        for i in m .. high:
          ei.v[i, j] += g * ei.ort[i]

proc cdiv[T](xr, xi, yr, yi: T): tuple[re, im: T] =
  # Complex scalar division.
  var r, d: T
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

proc hqr2[T](ei: var EigenvalueDecomposition[T]) =
  # Nonsymmetric reduction from Hessenberg to real Schur form.
  #
  # This is derived from the Algol procedure hqr2,
  # by Martin and Wilkinson, Handbook for Auto. Comp.,
  # Vol.ii-Linear Algebra, and the corresponding
  # Fortran subroutine in EISPACK.

  # Initialize
  let nn = ei.v.n
  var n = nn - 1
  let low = 0
  let high = nn - 1
  let eps = epsilon(T)
  var exshift = T(0)
  var p, q, r, s, z, t, w, x, y: T
  # Store roots isolated by balanc and compute matrix norm
  var norm = T(0)
  for i in 0 ..< nn:
    if i < low or i > high:
      ei.d[i] = ei.h[i, i]
      ei.e[i] = T(0)
    for j in max(i - 1, 0) ..< nn:
      norm = norm + abs(ei.h[i, j])
  # Outer loop over eigenvalue index
  var iter = 0
  while n >= low:
    # Look for single small sub-diagonal element
    var l = n
    while l > low:
      s = abs(ei.h[l - 1, l - 1]) + abs(ei.h[l, l])
      if s == 0:
        s = norm
      if abs(ei.h[l, l - 1]) < eps * s:
        break
      l.dec
    # Check for convergence
    # One root found
    if l == n:
      ei.h[n, n] = ei.h[n, n] + exshift
      ei.d[n] = ei.h[n, n]
      ei.e[n] = T(0)
      n.dec
      iter = 0
    # Two roots found
    elif l == n - 1:
      w = ei.h[n, n - 1] * ei.h[n - 1, n]
      p = (ei.h[n - 1, n - 1] - ei.h[n, n]) / T(2)
      q = p * p + w
      z = sqrt(abs(q))
      ei.h[n, n] = ei.h[n, n] + exshift
      ei.h[n - 1, n - 1] = ei.h[n - 1, n - 1] + exshift
      x = ei.h[n, n]
      # Real pair
      if q >= 0:
        if p >= 0:
          z = p + z
        else:
          z = p - z
        ei.d[n - 1] = x + z
        ei.d[n] = ei.d[n - 1]
        if z != 0:
          ei.d[n] = x - w / z
        ei.e[n - 1] = T(0)
        ei.e[n] = T(0)
        x = ei.h[n, n - 1]
        s = abs(x) + abs(z)
        p = x / s
        q = z / s
        r = sqrt(p * p + q * q)
        p = p / r
        q = q / r
        # Row modification
        for j in n - 1 ..< nn:
          z = ei.h[n - 1, j]
          ei.h[n - 1, j] = q * z + p * ei.h[n, j]
          ei.h[n, j] = q * ei.h[n, j] - p * z
        # Column modification
        for i in 0 .. n:
          z = ei.h[i, n - 1]
          ei.h[i, n - 1] = q * z + p * ei.h[i, n]
          ei.h[i, n] = q * ei.h[i, n] - p * z
        # Accumulate transformations
        for i in low .. high:
          z = ei.v[i, n - 1]
          ei.v[i, n - 1] = q * z + p * ei.v[i, n]
          ei.v[i, n] = q * ei.v[i, n] - p * z
      # Complex pair
      else:
        ei.d[n - 1] = x + p
        ei.d[n] = x + p
        ei.e[n - 1] = z
        ei.e[n] = -z
      n = n - 2
      iter = 0
    # No convergence yet
    else:
      # Form shift
      x = ei.h[n, n]
      y = T(0)
      w = T(0)
      if l < n:
        y = ei.h[n - 1, n - 1]
        w = ei.h[n, n - 1] * ei.h[n - 1, n]
      # Wilkinson's original ad hoc shift
      if iter == 10:
        exshift += x
        for i in low .. n:
          ei.h[i, i] -= x
        s = abs(ei.h[n, n - 1]) + abs(ei.h[n - 1, n - 2])
        y = T(0.75) * s
        x = y
        w = T(-0.4375) * s * s
      # MATLAB's new ad hoc shift
      if iter == 30:
        s = (y - x) / T(2)
        s = s * s + w
        if s > 0:
          s = sqrt(s)
          if y < x:
            s = -s
          s = x - w / ((y - x) / T(2) + s)
          for i in low .. n:
            ei.h[i, i] -= s
          exshift += s
          w = T(0.964)
          y = w
          x = y
      iter.inc # (Could check iteration count here.)
         # Look for two consecutive small sub-diagonal elements
      var m = n - 2
      while m >= l:
        z = ei.h[m, m]
        r = x - z
        s = y - z
        p = (r * s - w) / ei.h[m + 1, m] + ei.h[m, m + 1]
        q = ei.h[m + 1, m + 1] - z - r - s
        r = ei.h[m+2, m + 1]
        s = abs(p) + abs(q) + abs(r)
        p = p / s
        q = q / s
        r = r / s
        if m == l:
          break
        if abs(ei.h[m, m - 1]) * (abs(q) + abs(r)) <
              eps * (abs(p) * (abs(ei.h[m - 1, m - 1]) + abs(z) +
              abs(ei.h[m + 1, m + 1]))):
          break
        m.dec
      for i in m + 2 .. n:
        ei.h[i, i - 2] = T(0)
        if i > m + 2:
          ei.h[i, i - 3] = T(0)
      # Double QR step involving rows l:n and columns m:n
      for k in m .. n - 1:
        let notlast = k != n - 1
        if k != m:
          p = ei.h[k, k - 1]
          q = ei.h[k + 1, k - 1]
          r = if notlast: ei.h[k+2, k - 1] else: T(0)
          x = abs(p) + abs(q) + abs(r)
          if x == 0:
            continue
          p = p / x
          q = q / x
          r = r / x
        s = sqrt(p * p + q * q + r * r)
        if p < 0:
          s = -s
        if s != 0:
          if k != m:
            ei.h[k, k - 1] = -s * x
          elif l != m:
            ei.h[k, k - 1] = -ei.h[k, k - 1]
          p = p + s
          x = p / s
          y = q / s
          z = r / s
          q = q / p
          r = r / p
          # Row modification
          for j in k ..< nn:
            p = ei.h[k, j] + q * ei.h[k + 1, j]
            if notlast:
              p = p + r * ei.h[k + 2, j]
              ei.h[k + 2, j] = ei.h[k + 2, j] - p * z
            ei.h[k, j] = ei.h[k, j] - p * x
            ei.h[k + 1, j] = ei.h[k + 1, j] - p * y
          # Column modification
          for i in 0 .. min(n, k + 3):
            p = x * ei.h[i, k] + y * ei.h[i, k + 1]
            if notlast:
              p = p + z * ei.h[i, k + 2]
              ei.h[i, k+2] = ei.h[i, k + 2] - p * r
            ei.h[i, k] = ei.h[i, k] - p
            ei.h[i, k + 1] = ei.h[i, k + 1] - p * q
          # Accumulate transformations
          for i in low .. high:
            p = x * ei.v[i, k] + y * ei.v[i, k + 1]
            if notlast:
              p = p + z * ei.v[i, k + 2]
              ei.v[i, k + 2] = ei.v[i, k + 2] - p * r
            ei.v[i, k] = ei.v[i, k] - p
            ei.v[i, k + 1] = ei.v[i, k + 1] - p * q
        # (s != 0)
      # k loop
    # check convergence
  # while n >= low
  # Backsubstitute to find vectors of upper triangular form
  if norm == 0:
    return
  for n in countdown(nn - 1, 0):
    p = ei.d[n]
    q = ei.e[n]
    # Real vector
    if q == 0:
      var l = n
      ei.h[n, n] = T(1)
      for i in countdown(n - 1, 0):
        w = ei.h[i, i] - p
        r = T(0)
        for j in l .. n:
          r = r + ei.h[i, j] * ei.h[j, n]
        if ei.e[i] < 0:
          z = w
          s = r
        else:
          l = i
          if ei.e[i] == 0:
            if w != 0:
              ei.h[i, n] = -r / w
            else:
              ei.h[i, n] = -r / (eps * norm)
          # Solve real equations
          else:
            x = ei.h[i, i + 1]
            y = ei.h[i + 1, i]
            q = (ei.d[i] - p) * (ei.d[i] - p) + ei.e[i] * ei.e[i]
            t = (x * s - z * r) / q
            ei.h[i, n] = t
            if abs(x) > abs(z):
              ei.h[i + 1, n] = (-r - w * t) / x
            else:
              ei.h[i + 1, n] = (-s - y * t) / z
          # Overflow control
          t = abs(ei.h[i, n])
          if (eps * t) * t > 1:
            for j in i .. n:
              ei.h[j, n] = ei.h[j, n] / t
    # Complex vector
    elif q < 0:
      var l = n - 1
      # Last vector component imaginary so matrix is triangular
      if abs(ei.h[n, n - 1]) > abs(ei.h[n - 1, n]):
        ei.h[n - 1, n - 1] = q / ei.h[n, n - 1]
        ei.h[n - 1, n] = -(ei.h[n, n] - p) / ei.h[n, n - 1]
      else:
        let (cdivr, cdivi) = cdiv(T(0), -ei.h[n - 1, n], ei.h[n - 1, n - 1] - p, q)
        ei.h[n - 1, n - 1] = cdivr
        ei.h[n - 1, n] = cdivi
      ei.h[n, n - 1] = T(0)
      ei.h[n, n] = T(1)
      for i in countdown(n - 2, 0):
        var ra, sa, vr, vi: T
        ra = T(0)
        sa = T(0)
        for j in l .. n:
          ra = ra + ei.h[i, j] * ei.h[j, n - 1]
          sa = sa + ei.h[i, j] * ei.h[j, n]
        w = ei.h[i, i] - p
        if ei.e[i] < 0:
          z = w
          r = ra
          s = sa
        else:
          l = i
          if ei.e[i] == 0:
            let (cdivr, cdivi) = cdiv(-ra, -sa, w, q)
            ei.h[i, n - 1] = cdivr
            ei.h[i, n] = cdivi
          else:
            # Solve complex equations
            x = ei.h[i, i + 1]
            y = ei.h[i + 1, i]
            vr = (ei.d[i] - p) * (ei.d[i] - p) + ei.e[i] * ei.e[i] - q * q
            vi = (ei.d[i] - p) * T(2) * q
            if vr == 0 and vi == 0:
              vr = eps * norm * (abs(w) + abs(q) +
                   abs(x) + abs(y) + abs(z))
            let (cdivr, cdivi) = cdiv(x * r - z * ra + q * sa, x * s - z * sa -
                q * ra, vr, vi)
            ei.h[i, n - 1] = cdivr
            ei.h[i, n] = cdivi
            if abs(x) > (abs(z) + abs(q)):
              ei.h[i + 1, n - 1] = (-ra - w * ei.h[i, n - 1] + q * ei.h[i, n]) / x
              ei.h[i + 1, n] = (-sa - w * ei.h[i, n] - q * ei.h[i, n - 1]) / x
            else:
              let (cdivr, cdivi) = cdiv(-r - y * ei.h[i, n - 1], -s - y * ei.h[
                  i, n], z, q)
              ei.h[i + 1, n - 1] = cdivr
              ei.h[i + 1, n] = cdivi
          # Overflow control
          t = max(abs(ei.h[i, n - 1]), abs(ei.h[i, n]))
          if (eps * t) * t > 1:
            for j in i .. n:
              ei.h[j, n - 1] = ei.h[j, n - 1] / t
              ei.h[j, n] = ei.h[j, n] / t
  # Vectors of isolated roots
  for i in 0 ..< nn:
    if i < low or i > high:
      for j in i ..< nn:
        ei.v[i, j] = ei.h[i, j]
  # Back transformation to get eigenvectors of original matrix
  for j in countdown(nn - 1, low):
    for i in low .. high:
      z = T(0)
      for k in low .. min(j, high):
        z = z + ei.v[i, k] * ei.h[k, j]
      ei.v[i, j] = z

proc eig*[T](a: sink Matrix[T]): EigenvalueDecomposition[T] =
  ## Check for symmetry, then construct the eigenvalue decomposition.
  ##
  ## - parameter ``a``: Square matrix
  ## - ``return``: Structure to access D and V.
  # assert(a.n == a.m and a.n >= 1)
  let n = a.n
  result.d = newSeq[T](n)
  result.e = newSeq[T](n)
  var isSymmetric = true
  for j in 0 ..< n:
    if not isSymmetric:
      break
    for i in 0 ..< n:
      if not isSymmetric:
        break
      isSymmetric = a[i, j] == a[j, i]
  if isSymmetric:
    result.v = a
    # Tridiagonalize.
    result.tred2()
    # Diagonalize.
    result.tql2()
  else:
    result.h = a
    result.v = matrixUninit[T](n, n)
    result.ort = newSeq[T](n)
    # Reduce to Hessenberg form.
    result.orthes()
    # Reduce Hessenberg to real Schur form.
    result.hqr2()

proc getV*[T](ei: EigenvalueDecomposition[T]): lent Matrix[T] {.inline.} =
  ## Return the eigenvector matrix
  ei.v

proc getRealEigenvalues*[T](ei: EigenvalueDecomposition[T]): lent seq[T] {.inline.} =
  ## Return the real parts of the eigenvalues
  ##
  ## ``return``: ``real(diag(D))``
  ei.d

proc getImagEigenvalues*[T](ei: EigenvalueDecomposition[T]): lent seq[T] {.inline.} =
  ## Return the imaginary parts of the eigenvalues
  ##
  ## ``return``: ``imag(diag(D))``
  ei.e

proc getD*[T](ei: EigenvalueDecomposition[T]): Matrix[T] =
  ## Return the block diagonal eigenvalue matrix
  let n = ei.v.n
  result = matrixUninit[T](n, n)
  for i in 0 ..< n:
    for j in 0 ..< n:
      result[i, j] = T(0)
    result[i, i] = ei.d[i]
    if ei.e[i] > 0:
      result[i, i + 1] = ei.e[i]
    elif ei.e[i] < 0:
      result[i, i - 1] = ei.e[i]
