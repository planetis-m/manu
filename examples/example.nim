import math, times, strutils, manu

template ff(f: float, prec: int = 3): string =
  formatFloat(f, ffDecimal, prec)

proc magic(n: int): Matrix[float] =
  var m = matrix[float](n, n)
  # Odd order
  if n mod 2 == 1:
    let a = (n + 1) div 2
    let b = n + 1
    for j in 0 ..< n:
      for i in 0 ..< n:
        m[i, j] = float(n * ((i + j + a) mod n) + ((i + 2 * j + b) mod n) + 1)
  # Doubly Even Order
  elif n div 4 == 0:
    for j in 0 ..< n:
      for i in 0 ..< n:
        if (i + 1) div 2 mod 2 == (j + 1) div 2 mod 2:
          m[i, j] = float(n * n - n * i - j)
        else:
          m[i, j] = float(n * i + j + 1)
  # Singly Even Order
  else:
    let p = n div 2
    let k = (n - 2) div 4
    let a = magic(p)
    for j in 0 ..< p:
      for i in 0 ..< p:
        let aij = a[i, j]
        m[i, j] = aij
        m[i, j + p] = aij + float(2 * p * p)
        m[i + p, j] = aij + float(3 * p * p)
        m[i + p, j + p] = aij + float(p * p)
    for i in 0 ..< p:
      for j in 0 ..< k:
        swap(m[i, j], m[i + p, j])
      for j in n - k + 1 ..< n:
        swap(m[i, j], m[i + p, j])
    swap(m[k, 0], m[k + p, 0])
    swap(m[k, k], m[k + p, k])
  return m

proc main() =
  # Tests LU, QR, SVD and symmetric Eig decompositions.
  #
  #   n       = order of magic square.
  #   trace   = diagonal sum, should be the magic sum, (n^3 + n)/2.
  #   max_eig = maximum eigenvalue of (A + A')/2, should equal trace.
  #   rank    = linear algebraic rank,
  #             should equal n if n is odd, be less than n if n is even.
  #   cond    = L_2 condition number, ratio of singular values.
  #   lu_res  = test of LU factorization, norm1(L*U-A(p,:))/(n*eps).
  #   qr_res  = test of QR factorization, norm1(Q*R-A)/(n*eps).

  echo("    Test of Matrix object, using magic squares.")
  echo("    See MagicSquareExample.main() for an explanation.")
  echo("      n     trace       max_eig   rank        cond      lu_res      qr_res\n")

  let startTime = epochTime()
  let eps = pow(2.0, -52.0)
  var buf = ""
  for n in 3 .. 32:
    buf.setLen(0)
    buf.add(align($n, 7))

    let M = magic(n)
    let t = int(M.trace())
    buf.add(align($t, 10))

    let E = eig((M + M.transpose()) * 0.5)
    let d = E.getRealEigenvalues()
    buf.add(align(ff(d[n-1]), 14))

    let r = M.rank()
    buf.add(align($r, 7))

    let c = M.cond()
    if c < 1.0 / eps:
      buf.add(align(ff(c), 12))
    else:
      buf.add("         Inf")

    let LU = lu(M)
    let L = LU.getL()
    let U = LU.getU()
    let p = LU.getPivot()
    var R = L * U - M[p, 0 .. n-1]
    var res = R.norm1() / (float(n) * eps)
    buf.add(align(ff(res), 12))

    let QR = qr(M)
    let Q = QR.getQ()
    R = QR.getR()
    R = Q * R - M
    res = R.norm1() / (float(n) * eps)
    buf.add(align(ff(res), 12))

    echo buf

  let stopTime = epochTime()
  let etime = (stopTime - startTime) / 1000.0
  echo("\nElapsed Time = ", align(ff(etime), 12), " seconds")
  echo("Adios")

main()
