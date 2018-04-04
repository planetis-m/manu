import math

template checkBounds(cond: untyped, msg = "") =
   when compileOption("boundChecks"):
      {.line.}:
         if not cond:
            raise newException(IndexError, msg)

template newData() =
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

type Matrix* = object
   # Array for internal storage of elements.
   data: seq[seq[float]]
   # Row and column dimensions.
   m*, n*: int

proc matrix*(m, n: int): Matrix =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   newData()

proc matrix*(m, n: int, s: float): Matrix =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   newData()
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i][j] = s

proc matrix*(data: seq[seq[float]]): Matrix =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
   when compileOption("assertions"):
      for i in 0 ..< result.m:
         assert(data[i].len == result.n, "All rows must have the same length.")
   result.data = data

proc matrix*(data: seq[float], m: int): Matrix =
   ## Construct a matrix from a one-dimensional packed array.
   ## ``data`` is a one-dimensional array of float, packed by columns (ala Fortran).
   ## Array length must be a multiple of ``m``.
   result.m = m
   result.n = if m != 0: data.len div m else: 0
   assert result.m * result.n == data.len, "Array length must be a multiple of m."
   newData()
   for i in 0 ..< m:
      for j in 0 ..< result.n:
         result.data[i][j] = data[i + j * m]

proc getArray*(m: Matrix): seq[seq[float]] =
   ## Copy the internal two-dimensional array.
   result = m.data

proc getColumnPacked*(m: Matrix): seq[float] =
   ## Make a one-dimensional column packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i + j * m.m] = m.data[i][j]

proc getRowPacked*(m: Matrix): seq[float] =
   ## Make a one-dimensional row packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i * m.n + j] = m.data[i][j]

proc rowDimension*(m: Matrix): int =
   ## Get row dimension.
   m.m

proc columnDimension*(m: Matrix): int =
   ## Get column dimension.
   m.n

proc `[]`*(m: Matrix, i, j: int): float =
   ## Get a single element.
   m.data[i][j]

proc `[]`*(m: Matrix, r, c: Slice[int]): Matrix =
   ## Get a submatrix,
   ## ``m[i0 .. i1, j0 .. j1]``
   checkBounds(r.a >= 0 and r.b < m.m, "Submatrix dimensions")
   checkBounds(c.a >= 0 and c.b < m.n, "Submatrix dimensions")
   result.m = r.b - r.a + 1
   result.n = c.b - c.a + 1
   newData()
   for i in r.a .. r.b:
      for j in c.a .. c.b:
         result.data[i - r.a][j - c.a] = m.data[i][j]

proc `[]`*(m: Matrix, r, c: openarray[int]): Matrix =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result.m = r.len
   result.n = c.len
   newData()
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         result.data[i][j] = m.data[r[i]][c[j]]

proc `[]`*(m: Matrix, r: Slice[int], c: openarray[int]): Matrix =
   ## Get a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]]``
   checkBounds(r.a >= 0 and r.b < m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result.m = r.b - r.a + 1
   result.n = c.len
   newData()
   for i in r.a .. r.b:
      for j in 0 ..< c.len:
         result.data[i - r.a][j] = m.data[i][c[j]]

proc `[]`*(m: Matrix, r: openarray[int], c: Slice[int]): Matrix =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   checkBounds(c.a >= 0 and c.b < m.n, "Submatrix dimensions")
   result.m = r.len
   result.n = c.b - c.a + 1
   newData()
   for i in 0 ..< r.len:
      for j in c.a .. c.b:
         result.data[i][j - c.a] = m.data[r[i]][j]

proc `[]=`*(m: var Matrix, i, j: int, s: float) =
   ## Set a single element.
   m.data[i][j] = s

proc `[]=`*(m: var Matrix, r, c: Slice[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, j0 .. j1] = a``
   checkBounds(r.b - r.a + 1 == a.m, "Submatrix dimensions")
   checkBounds(c.b - c.a + 1 == a.n, "Submatrix dimensions")
   for i in r.a .. r.b:
      for j in c.a .. c.b:
         m.data[i][j] = a.data[i - r.a][j - c.a]

proc `[]=`*(m: var Matrix, r, c: openarray[int], a: Matrix) =
   ## Set a submatrix
   checkBounds(r.len == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         m.data[r[i]][c[j]] = a.data[i][j]

proc `[]=`*(m: var Matrix, r: openarray[int], c: Slice[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   checkBounds(c.b - c.a + 1 == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in c.a .. c.b:
         m.data[r[i]][j] = a.data[i][j - c.a]

proc `[]=`*(m: var Matrix, r: Slice[int], c: openarray[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]] = a``
   checkBounds(r.b - r.a + 1 == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in r.a .. r.b:
      for j in 0 ..< c.len:
         m.data[i][c[j]] = a.data[i - r.a][j]

proc `-`*(m: Matrix): Matrix =
   ## Unary minus
   result.m = m.m
   result.n = m.n
   newData()
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result.data[i][j] = -m.data[i][j]

proc `+`*(a, b: Matrix): Matrix =
   ## ``C = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result.m = a.m
   result.n = a.n
   newData()
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result.data[i][j] = a.data[i][j] + b.data[i][j]

proc `+=`*(a: var Matrix, b: Matrix) =
   ## ``A = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a.data[i][j] = a.data[i][j] + b.data[i][j]

proc `-`*(a, b: Matrix): Matrix =
   ## ``C = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result.m = a.m
   result.n = a.n
   newData()
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result.data[i][j] = a.data[i][j] + b.data[i][j]

proc `-=`*(a: var Matrix, b: Matrix) =
   ## ``A = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a.data[i][j] = a.data[i][j] + b.data[i][j]

proc `.*`*(a, b: Matrix): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result.m = a.m
   result.n = a.n
   newData()
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result.data[i][j] = a.data[i][j] * b.data[i][j]

proc `.*=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a.data[i][j] = a.data[i][j] * b.data[i][j]

proc `./`*(a, b: Matrix): Matrix =
   ## Element-by-element right division, ``C = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result.m = a.m
   result.n = a.n
   newData()
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result.data[i][j] = a.data[i][j] / b.data[i][j]

proc `./=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a.data[i][j] = a.data[i][j] / b.data[i][j]

proc `.\`*(a, b: Matrix): Matrix =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result.m = a.m
   result.n = a.n
   newData()
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         result.data[i][j] = b.data[i][j] / a.data[i][j]

proc `.\=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a.data[i][j] = b.data[i][j] / a.data[i][j]

proc `*`*(s: float, m: Matrix): Matrix =
   ## Multiply a matrix by a scalar, ``C = s*A``
   result.m = m.m
   result.n = m.n
   newData()
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result.data[i][j] = s * m.data[i][j]

proc `*=`*(m: var Matrix, s: float) =
   ## Multiply a matrix by a scalar in place, ``A = s*A``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m.data[i][j] = s * m.data[i][j]

proc `*`*(a, b: Matrix): Matrix =
   ## Linear algebraic matrix multiplication, ``A * B``
   assert(b.m == a.n, "Matrix inner dimensions must agree.")
   result.m = a.m
   result.n = b.n
   newData()
   var b_colj = newSeq[float](a.n)
   for j in 0 ..< b.n:
      for k in 0 ..< a.n:
         b_colj[k] = b.data[k][j]
      for i in 0 ..< a.m:
         var a_rowi = a.data[i]
         var s = 0.0
         for k in 0 ..< a.n:
            s += a_rowi[k] * b_colj[k]
         result.data[i][j] = s
