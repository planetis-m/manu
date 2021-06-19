## Singular Value Decomposition
## ============================
##
## For an m-by-n matrix A with m >= n, the singular value decomposition is
## an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
## an n-by-n orthogonal matrix V so that ``A = U*S*V'``.
##
## The singular values, ``sigma[k] = S[k, k]``, are ordered so that
## ``sigma[0] >= sigma[1] >= ... >= sigma[n-1]``.
##
## The singular value decompostion always exists, so the constructor will
## never fail.  The matrix condition number and the effective numerical
## rank can be computed from this decomposition.
import matrix, std/[math, fenv]

type
  SingularValueDecomposition*[T] = object
    u, v: Matrix[T] # Arrays for internal storage of U and V.
    s: seq[T]       # Array for internal storage of singular values.

proc svd*[T](a: sink Matrix[T]): SingularValueDecomposition[T] =
  ## Construct the singular value decomposition.
  ##
  ## ``return``: Structure to access U, S and V.
  # Derived from LINPACK code.
  # Initialize.
  let m = a.m
  let n = a.n
  #var a = a
  # assert(m >= n, "SVD only works for m >= n") and set nu = n
  let nu = min(m, n)
  result.s = newSeq[T](min(m + 1, n)) # n<=m so n is min
  result.u = matrix[T](m, nu)
  result.v = matrix[T](n, n)
  var e = newSeq[T](n)
  var work = newSeq[T](m)
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
      result.s[k] = T(0)
      for i in k ..< m:
        result.s[k] = hypot(result.s[k], a[i, k])
      if result.s[k] != 0:
        if a[k, k] < 0:
          result.s[k] = -result.s[k]
        for i in k ..< m:
          a[i, k] /= result.s[k]
        a[k, k] += T(1)
      result.s[k] = -result.s[k]
    for j in k + 1 ..< n:
      if k < nct and result.s[k] != 0:
        # Apply the transformation.
        var t = T(0)
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
      e[k] = T(0)
      for i in k + 1 ..< n:
        e[k] = hypot(e[k], e[i])
      if e[k] != 0:
        if e[k + 1] < 0:
          e[k] = -e[k]
        for i in k + 1 ..< n:
          e[i] /= e[k]
        e[k + 1] += T(1)
      e[k] = -e[k]
      if k + 1 < m and e[k] != 0:
        # Apply the transformation.
        for i in k + 1 ..< m:
          work[i] = T(0)
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
    result.s[nct] = a[nct, nct]
  if m < p:
    result.s[p - 1] = T(0)
  if nrt + 1 < p:
    e[nrt] = a[nrt, p - 1]
  e[p - 1] = T(0)
  # If required, generate U.
  if wantu:
    for j in nct ..< nu:
      for i in 0 ..< m:
        result.u[i, j] = T(0)
      result.u[j, j] = T(1)
    for k in countdown(nct - 1, 0):
      if result.s[k] != 0:
        for j in k + 1 ..< nu:
          var t = T(0)
          for i in k ..< m:
            t += result.u[i, k] * result.u[i, j]
          t = -t / result.u[k, k]
          for i in k ..< m:
            result.u[i, j] += t * result.u[i, k]
        for i in k ..< m:
          result.u[i, k] = -result.u[i, k]
        result.u[k, k] = T(1) + result.u[k, k]
        for i in 0 .. k - 2:
          result.u[i, k] = T(0)
      else:
        for i in 0 ..< m:
          result.u[i, k] = T(0)
        result.u[k, k] = T(1)
  # If required, generate V.
  if wantv:
    for k in countdown(n - 1, 0):
      if k < nrt and e[k] != 0:
        for j in k + 1 ..< nu:
          var t = T(0)
          for i in k + 1 ..< n:
            t += result.v[i, k] * result.v[i, j]
          t = -t / result.v[k + 1, k]
          for i in k + 1 ..< n:
            result.v[i, j] += t * result.v[i, k]
      for i in 0 ..< n:
        result.v[i, k] = T(0)
      result.v[k, k] = T(1)
  # Main iteration loop for the singular values.
  let pp = p - 1
  var iter = 0
  var eps = epsilon(T) # machine epsilon
  when T is float64:
    var tiny = pow(T(2), T(-966)) # no idea!
  elif T is float32:
    var tiny = pow(T(2), T(-120))
  while p > 0:
    var k = p - 2
    var kase = 0
    # Here is where a test for too many iterations would go.
    if iter == 500 or iter == 750:
      # echo("Svd taking a long time: making convergence criterion less exact.")
      eps = pow(T(0.8), eps)
      tiny = pow(T(0.8), tiny)
    assert(iter < 1000,
           "Svd not converging on matrix of size " & $m & " by " & $n)
    # This section of the program inspects for
    # negligible elements in the s and e arrays.  On
    # completion the variables kase and k are set as follows.

    # kase = 1     if s(p) and e[k-1] are negligible and k<p
    # kase = 2     if s(k) is negligible and k<p
    # kase = 3     if e[k-1] is negligible, k<p, and
    #              s(k), ..., s(p) are not negligible (qr step).
    # kase = 4     if e(p-1) is negligible (convergence).
    while k >= 0:
      if abs(e[k]) <=
            tiny + eps * (abs(result.s[k]) + abs(result.s[k + 1])):
        e[k] = T(0)
        break
      k.dec
    if k == p - 2:
      kase = 4
    else:
      var ks = p - 1
      while ks > k:
        var t = T(0)
        if ks != p: t += abs(e[ks])
        if ks != k + 1: t += abs(e[ks-1])
        if abs(result.s[ks]) <= tiny + eps * t:
          result.s[ks] = T(0)
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
      e[p - 2] = T(0)
      for j in countdown(p - 2, k):
        var t = hypot(result.s[j], f)
        let cs = result.s[j] / t
        let sn = f / t
        result.s[j] = t
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
      e[k - 1] = T(0)
      for j in k ..< p:
        var t = hypot(result.s[j], f)
        let cs = result.s[j] / t
        let sn = f / t
        result.s[j] = t
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
               abs(result.s[p - 1]), abs(result.s[p - 2])), abs(e[p - 2])),
               abs(result.s[k])), abs(e[k]))
      let sp = result.s[p - 1] / scale
      let spm1 = result.s[p - 2] / scale
      let epm1 = e[p - 2] / scale
      let sk = result.s[k] / scale
      let ek = e[k] / scale
      let b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / T(2)
      let c = (sp * epm1) * (sp * epm1)
      var shift = T(0)
      if b != 0 or c != 0:
        shift = sqrt(b * b + c)
        if b < 0:
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
          e[j - 1] = t
        f = cs * result.s[j] + sn * e[j]
        e[j] = cs * e[j] - sn * result.s[j]
        g = sn * result.s[j + 1]
        result.s[j + 1] = cs * result.s[j + 1]
        if wantv:
          for i in 0 ..< n:
            t = cs * result.v[i, j] + sn * result.v[i, j + 1]
            result.v[i, j + 1] = -sn * result.v[i, j] + cs * result.v[i, j + 1]
            result.v[i, j] = t
        t = hypot(f, g)
        cs = f / t
        sn = g / t
        result.s[j] = t
        f = cs * e[j] + sn * result.s[j + 1]
        result.s[j + 1] = -sn * e[j] + cs * result.s[j + 1]
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
      if result.s[k] <= 0:
        result.s[k] = if result.s[k] < 0: -result.s[k] else: T(0)
        if wantv:
          for i in 0 .. pp:
            result.v[i, k] = -result.v[i, k]
      # Order the singular values.
      while k < pp:
        if result.s[k] >= result.s[k + 1]:
          break
        swap(result.s[k + 1], result.s[k])
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

proc getU*[T](sv: SingularValueDecomposition[T]): lent Matrix[T] {.inline.} =
  ## Return the left singular vectors
  sv.u

proc getV*[T](sv: SingularValueDecomposition[T]): lent Matrix[T] {.inline.} =
  ## Return the right singular vectors
  sv.v

proc getSingularValues*[T](sv: SingularValueDecomposition[T]): lent seq[T] {.inline.} =
  ## Return the one-dimensional array of singular values.
  ##
  ## ``return``: diagonal of S
  sv.s

proc getS*[T](sv: SingularValueDecomposition[T]): Matrix[T] =
  ## Return the diagonal matrix of singular values.
  result = matrixUninit[T](sv.v.m, sv.v.n)
  for i in 0 ..< sv.v.m:
    for j in 0 ..< sv.v.n:
      result[i, j] = T(0)
    result[i, i] = sv.s[i]

proc norm2*[T](sv: SingularValueDecomposition[T]): T {.inline.} =
  ## Two norm.
  ##
  ## ``return``: max(S)
  sv.s[0]

proc cond*[T](sv: SingularValueDecomposition[T]): T {.inline.} =
  ## Two norm condition number.
  ##
  ## ``return``: max(S)/min(S)
  sv.s[0] / sv.s[^1]

proc rank*[T](sv: SingularValueDecomposition[T]): int =
  ## Effective numerical matrix rank.
  ##
  ## ``return``: Number of nonnegligible singular values.
  let eps = epsilon(T)
  let tol = T(max(sv.u.m, sv.v.m)) * sv.s[0] * eps
  for d in sv.s:
    if d > tol:
      result.inc
