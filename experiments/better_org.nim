import math

type
   Matrix = object
      # Array for internal storage of elements.
      data: seq[seq[float]]
      # Row and column dimensions.
      m, n: int
   CholeskyDecomposition = object
      # triangular factor
      l: Matrix
      # is symmetric and positive definite flag.
      isspd: bool

proc matrix(m, n: int): Matrix =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

proc matrix(data: seq[seq[float]]): Matrix =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
   result.data = data

proc `[]`(m: Matrix, i, j: int): float =
   ## Get a single element.
   m.data[i][j]

proc `[]`(m: var Matrix, i, j: int): var float =
   ## Get a single element.
   m.data[i][j]

proc `[]=`(m: var Matrix, i, j: int, s: float) =
   ## Set a single element.
   m.data[i][j] = s

proc mgetRow(m: var Matrix, i: int): var seq[float] =
   m.data[i]

proc chol(m: Matrix): CholeskyDecomposition =
   let n = m.m
   result.l = matrix(n, n)
   result.isspd = m.n == m.m
   # Main loop.
   for j in 0 ..< n:
      var lRowj = result.l.mgetRow(j)
      var d = 0.0
      for k in 0 ..< j:
         var lRowk = result.l.mgetRow(k)
         var s = 0.0
         for i in 0 ..< k:
            s += lRowk[i] * lRowj[i]
         s = (m[j, k] - s) / result.l[k, k]
         lRowj[k] = s
         d = d + s * s
         result.isspd = result.isspd and m[k, j] == m[j, k]
      d = m[j, j] - d
      result.isspd = result.isspd and d > 0.0
      result.l[j, j] = sqrt(max(d, 0.0))

let
   columnwise = @[@[18.0, 22, 54, 42], @[22.0, 70, 86, 62], @[54.0, 86, 174, 134], @[42.0, 62, 134, 106]]
   mat = matrix(columnwise)
   c = chol(mat)
