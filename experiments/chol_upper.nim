import math, random, times, strutils

type
   AnyMatrix {.explain.} = concept m, var mvar, type M
      M.m is int
      M.n is int
      m[int, int] is float
      mvar[int, int] = float

   MatrixA = object
      data: seq[seq[float]]
      m, n: int

   MatrixB = object
      data: seq[float]
      m, n: int

   CholeskyDecomposition[AnyMatrix] = object
      # triangular factor
      l: AnyMatrix
      # is symmetric and positive definite flag.
      isspd: bool

proc `[]`(m: MatrixA, i, j: int): float {.inline.} =
   m.data[i][j]

proc `[]`(m: var MatrixA, i, j: int): var float {.inline.} =
   m.data[i][j]

proc `[]=`(m: var MatrixA, i, j: int, v: float) {.inline.} =
   m.data[i][j] = v

proc rowAddr(m: var MatrixA, i: int): ptr seq[float] =
   m.data[i].addr

template newData() =
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

proc matrixA(m, n: int): MatrixA =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   newData()

proc randomMatrix(m, n: Natural): MatrixA =
   const maxVal = 1000
   result.m = m
   result.n = n
   newData()
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i][j] = rand(maxVal).float

proc getRowPacked(m: MatrixA): seq[float] =
   ## Make a one-dimensional row packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i * m.n + j] = m.data[i][j]

proc getColumnPacked(m: MatrixA): seq[float] =
   ## Make a one-dimensional column packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i + j * m.m] = m.data[i][j]

proc `[]`(m: MatrixB, i, j: int): float {.inline.} =
   m.data[i * m.n + j]

proc `[]`(m: var MatrixB, i, j: int): var float {.inline.} =
   m.data[i * m.n + j]

proc `[]=`(m: var MatrixB, i, j: int, v: float) {.inline.} =
   m.data[i * m.n + j] = v

proc matrixB(m, n: int): MatrixB =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   newSeq(result.data, m * n)

proc matrix(data: seq[float], m: int): MatrixB =
   result.m = m
   result.n = if m != 0: data.len div m else: 0
   result.data = data

proc cholL(a: AnyMatrix): CholeskyDecomposition[AnyMatrix] =
   ## Cholesky algorithm for symmetric and positive definite matrix.
   ## Structure to access L and isspd flag.
   ## parameter ``m``: Square, symmetric matrix.
   let n = a.m
   result.l = matrixA(n, n)
   result.isspd = a.n == a.m
   # Main loop.
   for j in 0 ..< n:
      var lRowj = result.l.rowAddr(j)
      var d = 0.0
      for k in 0 ..< j:
         var lRowk = result.l.rowAddr(k)
         var s = 0.0
         for i in 0 ..< k:
            s += lRowk[i] * lRowj[i]
         s = (a[j, k] - s) / result.l[k, k]
         lRowj[k] = s
         d = d + s * s
         result.isspd = result.isspd and a[k, j] == a[j, k]
      d = a[j, j] - d
      result.isspd = result.isspd and d > 0.0
      result.l[j, j] = sqrt(max(d, 0.0))

proc cholR(a: MatrixB): CholeskyDecomposition[MatrixB] =
   # Initialize.
   let n = a.m
   result.l = matrixB(n, n)
   result.isspd = a.n == a.m
   # Main loop.
   for j in 0 ..< n:
      var d = 0.0
      for k in 0 ..< j:
         var s = a[k, j]
         for i in 0 ..< k:
            s = s - result.l[i, k] * result.l[i, j]
         s = s / result.l[k, k]
         result.l[k, j] = s
         d = d + s * s
         result.isspd = result.isspd and a[k, j] == a[j, k] 
      d = a[j, j] - d
      result.isspd = result.isspd and d > 0.0
      result.l[j, j] = sqrt(max(d, 0.0))

proc main() =
   const n = 1000
   let a = randomMatrix(n, n)
   let b = matrix(a.getRowPacked(), n)
   block:
      # time Jama version
      let start = epochTime()
      let c = chol(a)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us naive storage"
   block:
      # time standard ijk
      let start = epochTime()
      let c = chol(b)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us packed storage"

main()
