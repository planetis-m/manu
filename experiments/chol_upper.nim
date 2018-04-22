import math, times, strutils
import matrix_versions

type
   CholeskyDecomposition[T] = object
      # triangular factor
      l: T
      # is symmetric and positive definite flag.
      isspd: bool

proc cholL(a: MatrixA): CholeskyDecomposition[MatrixA] =
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

proc cholL(a: MatrixB): CholeskyDecomposition[MatrixB] =
   ## Cholesky algorithm for symmetric and positive definite matrix.
   ## Structure to access L and isspd flag.
   ## parameter ``m``: Square, symmetric matrix.
   let n = a.m
   result.l = matrixB(n, n)
   result.isspd = a.n == a.m
   # Main loop.
   for j in 0 ..< n:
      var d = 0.0
      for k in 0 ..< j:
         var s = 0.0
         for i in 0 ..< k:
            s += result.l[i, i] * result.l[j, i]
         s = (a[j, k] - s) / result.l[k, k]
         result.l[j, k] = s
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
      let c = cholL(a)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us naive storage"
   block:
      # time standard ijk
      let start = epochTime()
      let c = cholR(b)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us packed storage"

main()
