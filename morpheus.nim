import math, random, strutils
import matrix, cholesky, qr, lu#, eigen, svd

export matrix, cholesky, qr, lu

proc t*(m: Matrix): Matrix =
   ## Matrix transpose
   result.m = m.n
   result.n = m.m
   newData()
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result.data[j][i] = m.data[i][j]

proc identity*(m, n: int): Matrix =
   ## Generate identity matrix,
   ## returns An m-by-n matrix with ones on the diagonal and zeros elsewhere.
   result.m = m
   result.n = n
   newData()
   for i in 0 ..< m:
      for j in 0 ..< n:
         if i == j:
            result.data[i][j] = 1.0

proc norm1*(m: Matrix): float =
   ## One norm,
   ## returns maximum column sum.
   for j in 0 ..< m.n:
      var s = 0.0
      for i in 0 ..< m.m:
         s += abs(m.data[i][j])
      result = max(result, s)

# proc norm2*(m: Matrix): float =
#    ## Two norm,
#    ## returns maximum singular value.
#    svd(m).norm2())

proc normInf*(m: Matrix): float =
   ## Infinity norm,
   ## returns maximum row sum.
   for i in 0 ..< m.m:
      var s = 0.0
      for j in 0 ..< m.n:
         s += abs(m.data[i][j])
      result = max(result, s)

proc normF*(m: Matrix): float =
   ## Frobenius norm,
   ## returns sqrt of sum of squares of all elements.
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result = hypot(result, m.data[i][j])

proc solve(a, b: Matrix): Matrix =
   ## Solve ``A*X = B``,
   ## returns solution if A is square, least squares solution otherwise
   ## parameter ``b``: the right hand side
   if m == n:
      lu(a).solve(b)
   else:
      qr(a).solve(b)

proc solveTranspose(a, b: Matrix): Matrix =
   ## Solve ``X*A = B``, which is also ``A'*X' = B'``,
   ## returns solution if A is square, least squares solution otherwise.
   ## parameter ``b``: the right hand side
   t(m).solve(b.t())

proc inverse(m: Matrix): Matrix =
   ## Matrix inverse or pseudoinverse,
   ## returns inverse(A) if A is square, pseudoinverse otherwise.
   solve(m, identity(m.m, m.m))

proc det(m: Matrix): float =
   ## Matrix determinant
   lu(m).det()

# proc rank(m: Matrix): int =
#    ## Matrix rank,
#    ## returns effective numerical rank, obtained from SVD.
#    svd(m).rank()

# proc cond(m: Matrix): float =
#    ## Matrix condition (2 norm),
#    ## returns ratio of largest to smallest singular value.
#    svd(m).cond()

proc trace*(m: Matrix): float =
   ## Matrix trace,
   ## returns the sum of the diagonal elements.
   for i in 0 ..< min(m.m, m.n):
      result += m.data[i][i]

proc randMatrix*(m, n: int): Matrix =
   ## Generate matrix with random elements,
   ## returns an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   newData()
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i][j] = rand(1.0)

proc columnFormat(s: seq[float]): seq[string] =
   result = newSeq[string](s.len)
   for i, v in s:
      result[i] = formatFloat(v, ffDecimal, 6)
   var lenLeft = newSeq[int](s.len)
   var maxLenLeft = 0
   for i, f in result:
      let index = f.find('.')
      lenLeft[i]  = index
      maxLenLeft = max(maxLenLeft, lenLeft[i])
   for i in 0 ..< s.len:
      result[i] = spaces(maxLenLeft  - lenLeft[i]) & result[i]

proc `$`*(m: Matrix): string =
   var cols: seq[seq[string]]
   newSeq(cols, m.m)
   for i in 0 ..< m.m:
      cols[i] = columnFormat(m.data[i])
   result = ""
   for j in 0 ..< m.n:
      if j == 0:
         result.add "⎡"
      elif j == m.n - 1:
         result.add "⎣"
      else:
         result.add "⎢"
      for i in 0 ..< m.m:
         if i != 0:
            result.add "  "
         result.add cols[i][j]
      if j == 0:
         result.add "⎤\n"
      elif j == m.n - 1:
         result.add "⎦\n"
      else:
         result.add "⎥\n"
