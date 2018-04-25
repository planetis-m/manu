import times, strutils
import matrix_versions

proc multiply(a, b: MatrixA): MatrixA =
   result = matrixA(a.m, b.n)
   var bColj = newSeq[float](a.n)
   for j in 0 ..< b.n:
      for k in 0 ..< a.n:
         bColj[k] = b[k, j]
      for i in 0 ..< a.m:
         var aRowi = a.rowUnsafeAddr(i)
         var s = 0.0
         for k in 0 ..< a.n:
            s += aRowi[k] * bColj[k]
         result[i, j] = s

proc multiply(a, b: MatrixB): MatrixB =
   result = matrixB(a.m, b.n)
   var bColj = newSeq[float](a.n)
   for i in 0 ..< b.n:
      for k in 0 ..< a.n:
         bColj[k] = b[k, i]
      for j in 0 ..< a.m:
         var s = 0.0
         for k in 0 ..< a.n:
            s += a[i, k] * bColj[k]
         result[i, j] = s

proc multiply(a, b: MatrixC): MatrixC =
   result = matrixC(a.m, b.n)
   var aRowi = newSeq[float](a.n)
   for j in 0 ..< a.m:
      for k in 0 ..< a.n:
         aRowi[k] = a[j, k]
      for i in 0 ..< b.n:
         var s = 0.0
         for k in 0 ..< a.n:
            s += aRowi[k] * b[k, j]
         result[i, j] = s

proc main() =
   const n = 1000
   let a = randomMatrix(n, n)
   let b = randomMatrix(n, n)
   let d = matrixB(a.getColumnPacked(), n)
   let e = matrixB(b.getColumnPacked(), n)
   let f = matrixC(a.getRowPacked(), n)
   let g = matrixC(b.getRowPacked(), n)
   block:
      # time Jama version
      let start = epochTime()
      let c =  multiply(a, b)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us naive storage"
   block:
      # time standard ijk
      let start = epochTime()
      let c = multiply(d, e)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us packed storage rowwise"
   block:
      # time standard ijk
      let start = epochTime()
      let c = multiply(f, g)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us packed storage columnwise"

main()
