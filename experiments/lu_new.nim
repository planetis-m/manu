import random, times, strutils

type
   MatrixA = object
      data: seq[seq[float]]
      m, n: int

   MatrixB = object
      data: seq[float]
      m, n: int

   LUDecomposition[T] = object
      # Array for internal storage of decomposition.
      lu: T
      # Internal storage of pivot vector.
      piv: seq[int]
      # Pivot sign.
      pivsign: int

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

proc `[]`(m: MatrixB, i, j: int): float {.inline.} =
   m.data[i * m.n + j]

proc `[]`(m: var MatrixB, i, j: int): var float {.inline.} =
   m.data[i * m.n + j]

proc `[]=`(m: var MatrixB, i, j: int, v: float) {.inline.} =
   m.data[i * m.n + j] = v

proc matrix(data: seq[float], m: int): MatrixB =
   result.m = m
   result.n = if m != 0: data.len div m else: 0
   result.data = data

proc lu(a: MatrixA): LUDecomposition[MatrixA] =
   ## LU Decomposition
   ## Structure to access L, U and piv.
   ## param ``a``: Rectangular matrix
   # Use a "left-looking", dot-product, Crout/Doolittle algorithm.
   let m = a.m
   let n = a.n
   result.lu = a
   result.piv = newSeq[int](m)
   for i in 0 ..< m:
      result.piv[i] = i
   result.pivsign = 1
   var luColj = newSeq[float](m)
   # Outer loop.
   for j in 0 ..< n:
      # Make a copy of the j-th column to localize references.
      for i in 0 ..< m:
         luColj[i] = result.lu[i, j]
      # Apply previous transformations.
      for i in 0 ..< m:
         var luRowi = result.lu.rowAddr(i)
         # Most of the time is spent in the following dot product.
         let kmax = min(i, j)
         var s = 0.0
         for k in 0 ..< kmax:
            s += luRowi[k] * luColj[k]
         luColj[i] -= s
         luRowi[j] = luColj[i]
      # Find pivot and exchange if necessary.
      var p = j
      for i in j + 1 ..< m:
         if abs(luColj[i]) > abs(luColj[p]):
            p = i
      if p != j:
         for k in 0 ..< n:
            swap(result.lu[p, k], result.lu[j, k])
         swap(result.piv[p], result.piv[j])
         result.pivsign = -result.pivsign
      # Compute multipliers.
      if j < m and result.lu[j, j] != 0.0:
         for i in j + 1 ..< m:
            result.lu[i, j] /= result.lu[j, j]

proc lu(a: MatrixB): LUDecomposition[MatrixB] =
   ## LU Decomposition, computed by Gaussian elimination.
   ##
   ## This constructor computes L and U with the "daxpy"-based elimination
   ## algorithm used in LINPACK and MATLAB.  In Java, we suspect the dot-product,
   ## Crout algorithm will be faster.  We have temporarily included this
   ## constructor until timing experiments confirm this suspicion.
   ##
   ## parameter ``A``: Rectangular matrix
   ## return: Structure to access L, U and piv.
   # Initialize.
   let m = a.m
   let n = a.n
   result.lu = a
   result.piv = newSeq[int](m)
   for i in 0 ..< m:
      result.piv[i] = i
   result.pivsign = 1
   # Main loop.
   for k in 0 ..< n:
      # Find pivot.
      var p = k
      for i in k + 1 ..< m:
         if abs(result.lu[i, k]) > abs(result.lu[p, k]):
            p = i
      # Exchange if necessary.
      if p != k:
         for j in 0 ..< n:
            swap(result.lu[p, j], result.lu[k, j])
         swap(result.piv[p], result.piv[k])
         result.pivsign = -result.pivsign
      # Compute multipliers and eliminate k-th column.
      if result.lu[k, k] != 0.0:
         for i in k + 1 ..< m:
            result.lu[i, k] /= result.lu[k, k]
            for j in k + 1 ..< n:
               result.lu[i, j] -= result.lu[i, k] * result.lu[k, j]

proc main() =
   const n = 1000
   let a = randomMatrix(n, n)
   let b = matrix(a.getRowPacked(), n)
   block:
      # time Jama version
      let start = epochTime()
      let c = lu(a)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us naive storage"
   block:
      # time standard ijk
      let start = epochTime()
      let c = lu(b)
      let duration = epochTime() - start
      echo formatFloat(duration, ffDecimal, 3), " us packed storage"

main()
