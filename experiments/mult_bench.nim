import random, times, strutils

type
   MatrixA = object
      data: seq[seq[float]]
      m, n: int

   MatrixB = object
      data: seq[float]
      m, n: int

proc `[]`(m: MatrixB, i, j: int): float {.inline.} =
   m.data[i * m.n + j]

proc `[]=`(m: var MatrixB, i, j: int, v: float) {.inline.} =
   m.data[i * m.n + j] = v

proc matrix(data: seq[float], m: int): MatrixB =
   result.m = m
   result.n = if m != 0: data.len div m else: 0
   result.data = data

proc multiply(a, b: MatrixB): MatrixB =
   result.m = a.m
   result.n = b.n
   newSeq(result.data, result.m * result.n)
   for i in 0 ..< b.n:
      for j in 0 ..< a.m:
         var s = 0.0
         for k in 0 ..< a.n:
            s += a[i, k] * b[k, j]
         result[i, j] = s

template newData() =
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

proc randomMatrix(m, n: Natural): MatrixA =
   const maxVal = 1000
   result.m = m
   result.n = n
   newData()
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i][j] = rand(maxVal).float

proc multiply(a, b: MatrixA): MatrixA =
   result.m = a.m
   result.n = b.n
   newData()
   var b_colj = newSeq[float](a.n)
   for j in 0 ..< b.n:
      for k in 0 ..< a.n:
         b_colj[k] = b.data[k][j]
      for i in 0 ..< a.m:
         var a_rowi = unsafeAddr a.data[i]
         var s = 0.0
         for k in 0 ..< a.n:
            s += a_rowi[k] * b_colj[k]
         result.data[i][j] = s

proc getRowPacked(m: MatrixA): seq[float] =
   ## Make a one-dimensional row packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i * m.n + j] = m.data[i][j]

proc main() =
   const n = 1000
   let a = randomMatrix(n, n)
   let b = randomMatrix(n, n)
   let d = matrix(a.getRowPacked(), n)
   let e = matrix(b.getRowPacked(), n)
   block:
      # time Jama version
      let start = epochTime()
      let c =  multiply(a, b)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us naive storage - optimized ijk algorithm"
   block:
      # time standard ijk
      let start = epochTime()
      let f = multiply(d, e)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us packed storage - ijk algorithm"

main()
