import math, morpheus/matrix

template newData() =
   newSeq(result.data, result.n)
   for i in 0 ..< result.n:
      newSeq(result.data[i], result.n)

type CholeskyDecomposition = object
   data: seq[seq[float]]
   n: int
   isspd: bool

proc cholL(m: Matrix): CholeskyDecomposition =
   result.n = m.m
   newData()
   result.isspd = m.n == m.m
   # Main loop.
   for j in 0 ..< result.n:
      var lRowj = addr result.data[j]
      var d = 0.0
      for k in 0 ..< j:
         var lRowk = addr result.data[k]
         var s = 0.0
         for i in 0 ..< k:
            s = s + lRowk[i] * lRowj[i]
         s = (m.data[j][k] - s) / result.data[k][k]
         lRowj[k] = s
         d = d + s * s
         result.isspd = result.isspd and m.data[k][j] == m.data[j][k]

      d = m.data[j][j] - d
      result.isspd = result.isspd and d > 0.0
      result.data[j][j] = sqrt(max(d, 0.0))

proc cholR(m: Matrix): CholeskyDecomposition =
   # Initialize.
   result.n = m.m
   newData()
   result.isspd = m.n == m.m
   # Main loop.
   for j in 0 ..< result.n:
      var d = 0.0
      for k in 0 ..< j:
         var s = m.data[k][j]
         for i in 0 ..< k:
            s = s - result.data[i][k] * result.data[i][j]
         s = s / result.data[k][k]
         result.data[k][j] = s
         d = d + s * s
         result.isspd = result.isspd and m.data[k][j] == m.data[j][k] 

      d = m.data[j][j] - d
      result.isspd = result.isspd and d > 0.0
      result.data[j][j] = sqrt(max(d, 0.0))

proc getFactor(c: CholeskyDecomposition): Matrix =
   result.data = c.data
   result.m = c.n
   result.n = c.n

proc check(a, b: Matrix) =
   let eps = pow(2.0, -52.0)
   let x_norm1 = a.norm1()
   let y_norm1 = b.norm1()
   let xmiy_norm1 = norm1(a - b)
   if x_norm1 == 0.0 and y_norm1 < 10.0 * eps: return
   if y_norm1 == 0.0 and x_norm1 < 10.0 * eps: return
   if xmiy_norm1 > 1000.0 * eps * max(x_norm1, y_norm1):
      raise newException(AssertionError, "The norm of (a-b) is too large: " & $xmiy_norm1)

let
   columnwise = @[25.0, 15, -5, 15, 18, 0, -5, 0, 11]
   mat = matrix(columnwise, 3)
   l = cholL(mat).getFactor
   r = cholR(mat).getFactor

check(l * l.transpose, r.transpose * r)
