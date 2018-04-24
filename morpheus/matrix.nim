import math, random, strutils

template checkBounds(cond: untyped, msg = "") =
   when compileOption("boundChecks"):
      {.line.}:
         if not cond:
            raise newException(IndexError, msg)

type Matrix* = object
   # Array for internal storage of elements.
   data: seq[float]
   # Row and column dimensions.
   m*, n*: int

proc matrix*(m, n: int): Matrix {.inline.} =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   result.data = newSeqUninitialized[float](m * n)

proc matrix*(m, n: int, s: float): Matrix =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   result.data = newSeqUninitialized[float](m * n)
   for i in 0 ..< result.data.len:
      result.data[i] = s

proc matrix*(data: seq[seq[float]]): Matrix =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
   when compileOption("assertions"):
      for i in 0 ..< result.m:
         assert(data[i].len == result.n, "All rows must have the same length.")
   result.data = newSeqUninitialized[float](result.m * result.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result.data[i * result.n + j] = data[i][j]

proc matrix*(data: seq[seq[float]], m, n: int): Matrix =
   ## Construct a matrix quickly without checking arguments.
   result.m = m
   result.n = n
   result.data = newSeqUninitialized[float](m * n)
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i][j]

proc matrix*(data: seq[float], m: int): Matrix =
   ## Construct a matrix from a one-dimensional packed array.
   ##
   ## parameter ``data``: one-dimensional array of float, packed by columns (ala Fortran).
   ## Array length must be a multiple of ``m``.
   let n = if m != 0: data.len div m else: 0
   result.m = m
   result.n = n
   assert(result.m * result.n == data.len, "Array length must be a multiple of m.")
   result.data = newSeqUninitialized[float](data.len)
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i + j * m]

proc randMatrix*(m, n: int): Matrix =
   ## Generate matrix with random elements.
   ##
   ## ``return``: an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   result.data = newSeqUninitialized[float](result.m * result.n)
   for i in 0 ..< result.data.len:
      result.data[i] = rand(1.0)

proc getArray*(m: Matrix): seq[seq[float]] =
   ## Make a two-dimensional array copy of the internal array.
   result = newSeq[seq[float]](m.m)
   for i in 0 ..< m.m:
      result[i] = newSeq[float](m.n)
      for j in 0 ..< m.n:
         result[i][j] = m.data[i * m.n + j]

proc getColumnPacked*(m: Matrix): seq[float] =
   ## Make a one-dimensional column packed copy of the internal array.
   result = newSeq[float](m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i + j * m.m] = m.data[i * m.n + j]

proc getRowPacked*(m: Matrix): seq[float] {.inline.} =
   ## Copy the internal one-dimensional row packed array.
   m.data

proc rowDimension*(m: Matrix): int {.inline.} =
   ## Get row dimension.
   m.m

proc columnDimension*(m: Matrix): int {.inline.} =
   ## Get column dimension.
   m.n

proc `[]`*(m: Matrix, i, j: int): float {.inline.} =
   ## Get a single element.
   m.data[i * m.n + j]

proc `[]`*(m: var Matrix, i, j: int): var float {.inline.} =
   ## Get a single element.
   m.data[i * m.n + j]

proc `[]=`*(m: var Matrix, i, j: int, s: float) =
   ## Set a single element.
   m.data[i * m.n + j] = s

proc `[]`*(m: Matrix, r, c: Slice[int]): Matrix =
   ## Get a submatrix,
   ## ``m[i0 .. i1, j0 .. j1]``
   checkBounds(r.a >= 0 and r.b < m.m, "Submatrix dimensions")
   checkBounds(c.a >= 0 and c.b < m.n, "Submatrix dimensions")
   result = matrix(r.b - r.a + 1, c.b - c.a + 1)
   for i in r.a .. r.b:
      for j in c.a .. c.b:
         result[i - r.a, j - c.a] = m[i, j]

proc `[]`*(m: Matrix, r, c: openarray[int]): Matrix =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result = matrix(r.len, c.len)
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         result[i, j] = m[r[i], c[j]]

proc `[]`*(m: Matrix, r: Slice[int], c: openarray[int]): Matrix =
   ## Get a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]]``
   checkBounds(r.a >= 0 and r.b < m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result = matrix(r.b - r.a + 1, c.len)
   for i in r.a .. r.b:
      for j in 0 ..< c.len:
         result[i - r.a, j] = m[i, c[j]]

proc `[]`*(m: Matrix, r: openarray[int], c: Slice[int]): Matrix =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   checkBounds(c.a >= 0 and c.b < m.n, "Submatrix dimensions")
   result = matrix(r.len, c.b - c.a + 1)
   for i in 0 ..< r.len:
      for j in c.a .. c.b:
         result[i, j - c.a] = m[r[i], j]

proc `[]=`*(m: var Matrix, r, c: Slice[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, j0 .. j1] = a``
   checkBounds(r.b - r.a + 1 == a.m, "Submatrix dimensions")
   checkBounds(c.b - c.a + 1 == a.n, "Submatrix dimensions")
   for i in r.a .. r.b:
      for j in c.a .. c.b:
         m[i, j] = a[i - r.a, j - c.a]

proc `[]=`*(m: var Matrix, r, c: openarray[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         m[r[i], c[j]] = a[i, j]

proc `[]=`*(m: var Matrix, r: openarray[int], c: Slice[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   checkBounds(c.b - c.a + 1 == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in c.a .. c.b:
         m[r[i], j] = a[i, j - c.a]

proc `[]=`*(m: var Matrix, r: Slice[int], c: openarray[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]] = a``
   checkBounds(r.b - r.a + 1 == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in r.a .. r.b:
      for j in 0 ..< c.len:
         m[i, c[j]] = a[i - r.a, j]

proc `-`*(m: Matrix): Matrix =
   ## Unary minus
   result = matrix(m.m, m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i, j] = -m[i, j]

proc `+`*(a, b: Matrix): Matrix =
   ## ``C = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = matrix(a.m, a.n)
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result[i, j] = a[i, j] + b[i, j]

proc `+=`*(a: var Matrix, b: Matrix) =
   ## ``A = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[i, j]

proc `-`*(a, b: Matrix): Matrix =
   ## ``C = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = matrix(a.m, a.n)
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result[i, j] = a[i, j] - b[i, j]

proc `-=`*(a: var Matrix, b: Matrix) =
   ## ``A = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[i, j]

proc `.*`*(a, b: Matrix): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = matrix(a.m, a.n)
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result[i, j] = a[i, j] * b[i, j]

proc `.*=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[i, j]

proc `./`*(a, b: Matrix): Matrix =
   ## Element-by-element right division, ``C = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = matrix(a.m, a.n)
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result[i, j] = a[i, j] / b[i, j]

proc `./=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[i, j]

proc `.\`*(a, b: Matrix): Matrix =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = matrix(a.m, a.n)
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result[i, j] = b[i, j] / a[i, j]

proc `.\=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[i, j] / a[i, j]

proc `*`*(m: Matrix, s: float): Matrix =
   ## Multiply a matrix by a scalar, ``C = s*A``
   result = matrix(m.m, m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i, j] = s * m[i, j]

proc `*=`*(m: var Matrix, s: float) =
   ## Multiply a matrix by a scalar in place, ``A = s*A``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = s * m[i, j]

proc `*`*(a, b: Matrix): Matrix =
   ## Linear algebraic matrix multiplication, ``A * B``
   assert(b.m == a.n, "Matrix inner dimensions must agree.")
   result = matrix(a.m, b.n)
   var bColj = newSeq[float](b.m)
   for j in 0 ..< b.n:
      for k in 0 ..< a.n:
         bColj[k] = b[k, j]
      for i in 0 ..< a.m:
         var s = 0.0
         for k in 0 ..< a.n:
            s += a[i, k] * bColj[k]
         result[i, j] = s

proc transpose*(m: Matrix): Matrix =
   ## Matrix transpose
   result = matrix(m.n, m.m)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[j, i] = m[i, j]

proc identity*(m, n: int): Matrix =
   ## Generate identity matrix.
   ##
   ## ``return``: An m-by-n matrix with ones on the diagonal and zeros elsewhere.
   result = matrix(m, n)
   for i in 0 ..< m:
      for j in 0 ..< n:
         if i == j:
            result[i, j] = 1.0
         else:
            result[i, j] = 0.0

proc norm1*(m: Matrix): float =
   ## One norm.
   ##
   ## ``return``: maximum column sum
   for j in 0 ..< m.n:
      var s = 0.0
      for i in 0 ..< m.m:
         s += abs(m[i, j])
      result = max(result, s)

proc normInf*(m: Matrix): float =
   ## Infinity norm.
   ##
   ## ``return``: maximum row sum
   for i in 0 ..< m.m:
      var s = 0.0
      for j in 0 ..< m.n:
         s += abs(m[i, j])
      result = max(result, s)

proc normF*(m: Matrix): float =
   ## Frobenius norm.
   ##
   ## ``return``: sqrt of sum of squares of all elements.
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result = hypot(result, m[i, j])

proc trace*(m: Matrix): float =
   ## Matrix trace.
   ##
   ## ``return``: the sum of the diagonal elements
   for i in 0 ..< min(m.m, m.n):
      result += m[i, i]

proc columnFormat(s: seq[float]): seq[string] =
   result = newSeq[string](s.len)
   var maxLen = 0
   for i, f in s:
      result[i] = formatEng(f)
      maxLen = max(maxLen, result[i].len)
   for i, f in result.mpairs:
      f = spaces(maxLen - f.len) & f

proc `$`*(m: Matrix): string =
   var cols: seq[seq[string]]
   newSeq(cols, m.n)
   var mColj = newSeq[float](m.m)
   for j in 0 ..< m.n:
      for i in 0 ..< m.m:
         mColj[i] = m[i, j]
      cols[j] = columnFormat(mColj)
   result = ""
   for i in 0 ..< m.m:
      if i == 0:
         result.add "⎡"
      elif i == m.m - 1:
         result.add "⎣"
      else:
         result.add "⎢"
      for j in 0 ..< m.n:
         if j != 0:
            result.add "  "
         result.add cols[j][i]
      if i == 0:
         result.add "⎤\n"
      elif i == m.m - 1:
         result.add "⎦"
      else:
         result.add "⎥\n"
