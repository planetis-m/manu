import math, random, strutils

template checkBounds(cond: untyped, msg = "") =
   when compileOption("boundChecks"):
      {.line.}:
         if not cond:
            raise newException(IndexError, msg)

type
   Matrix* = object
      m, n: int # Row and column dimensions.
      data: seq[float] # Array for internal storage of elements.

proc matrix*(m, n: int): Matrix {.inline.} =
   ## Construct an m-by-n matrix of zeros.
   result.m = m
   result.n = n
   result.data = newSeq[float](m * n)

proc matrix*(m, n: int, s: float): Matrix =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   result.data = newSeqUninitialized[float](m * n)
   for i in 0 ..< result.data.len:
      result.data[i] = s

template ones*(m, n: int): Matrix = matrix(m, n, 1.0)
template zeros*(m, n: int): Matrix = matrix(m, n)

proc matrix*(data: seq[seq[float]]): Matrix =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
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
   assert(m * n == data.len, "Array length must be a multiple of m.")
   result.m = m
   result.n = n
   result.data = newSeqUninitialized[float](data.len)
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i + j * m]

proc matrix*(n: int, data: sink seq[float]): Matrix =
   ## Construct a matrix from a one-dimensional packed array.
   ##
   ## parameter ``data``: one-dimensional array of float, packed by rows.
   ## Array length must be a multiple of ``n``.
   let m = if n != 0: data.len div n else: 0
   assert(m * n == data.len, "Array length must be a multiple of m.")
   result.m = m
   result.n = n
   result.data = data

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

proc dim*(m: Matrix): (int, int) {.inline.} =
   ## Get (row, column) dimensions tuple.
   (m.m, m.n)

proc m*(m: Matrix): int {.inline.} =
   ## Get row dimension.
   m.m

proc n*(m: Matrix): int {.inline.} =
   ## Get column dimension.
   m.n

proc rowDimension*(m: Matrix): int {.inline.} = m.m
proc columnDimension*(m: Matrix): int {.inline.} = m.n

proc `[]`*(m: Matrix, i, j: int): float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   m.data[i * m.n + j]

proc `[]`*(m: var Matrix, i, j: int): var float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   m.data[i * m.n + j]

proc `[]=`*(m: var Matrix, i, j: int, s: float) {.inline.} =
   ## Set a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   m.data[i * m.n + j] = s

template `^^`(dim, i: untyped): untyped =
  (when i is BackwardsIndex: dim - int(i) else: int(i))

proc `[]`*[U, V, W, X](m: Matrix, r: HSlice[U, V], c: HSlice[W, X]): Matrix =
   ## Get a submatrix,
   ## ``m[i0 .. i1, j0 .. j1]``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(ra >= 0 and rb < m.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(ca >= 0 and cb < m.n, "Submatrix dimensions")
   result = matrix(rb - ra + 1, cb - ca + 1)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = m[i + ra, j + ca]

proc `[]`*(m: Matrix, r, c: openarray[int]): Matrix =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result = matrix(r.len, c.len)
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         result[i, j] = m[r[i], c[j]]

proc `[]`*[U, V](m: Matrix, r: HSlice[U, V], c: openarray[int]): Matrix =
   ## Get a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]]``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(ra >= 0 and rb < m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result = matrix(rb - ra + 1, c.len)
   for i in 0 ..< result.m:
      for j in 0 ..< c.len:
         result[i, j] = m[i + ra, c[j]]

proc `[]`*[U, V](m: Matrix, r: openarray[int], c: HSlice[U, V]): Matrix =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(ca >= 0 and cb < m.n, "Submatrix dimensions")
   result = matrix(r.len, cb - ca + 1)
   for i in 0 ..< r.len:
      for j in 0 ..< result.n:
         result[i, j] = m[r[i], j + ca]

proc `[]=`*[U, V, W, X](m: var Matrix, r: HSlice[U, V], c: HSlice[W, X], a: Matrix) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, j0 .. j1] = a``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(rb - ra + 1 == a.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(cb - ca + 1 == a.n, "Submatrix dimensions")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         m[i + ra, j + ca] = a[i, j]

proc `[]=`*(m: var Matrix, r, c: openarray[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         m[r[i], c[j]] = a[i, j]

proc `[]=`*[U, V](m: var Matrix, r: HSlice[U, V], c: openarray[int], a: Matrix) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]] = a``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(rb - ra + 1 == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in 0 ..< a.m:
      for j in 0 ..< c.len:
         m[i + ra, c[j]] = a[i, j]

proc `[]=`*[U, V](m: var Matrix, r: openarray[int], c: HSlice[U, V], a: Matrix) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(cb - ca + 1 == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in 0 ..< a.n:
         m[r[i], j + ca] = a[i, j]

proc `-`*(m: sink Matrix): Matrix =
   ## Unary minus
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = -result[i, j]

proc `+`*(a: sink Matrix, b: Matrix): Matrix =
   ## ``C = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[i, j]

proc `+=`*(a: var Matrix, b: Matrix) =
   ## ``A = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[i, j]

proc `-`*(a: sink Matrix, b: Matrix): Matrix =
   ## ``C = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[i, j]

proc `-=`*(a: var Matrix, b: Matrix) =
   ## ``A = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[i, j]

proc `.*`*(a: sink Matrix, b: Matrix): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[i, j]

proc `.*=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[i, j]

proc `./`*(a: sink Matrix, b: Matrix): Matrix =
   ## Element-by-element right division, ``C = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / b[i, j]

proc `./=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[i, j]

proc `.\`*(a: sink Matrix, b: Matrix): Matrix =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[i, j] / result[i, j]

proc `.\=`*(a: var Matrix, b: Matrix) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[i, j] / a[i, j]

proc `+`*(m: sink Matrix, s: float): Matrix =
   ## Add a matrix to a scalar, ``C = s+A``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = s + result[i, j]

template `+`*(s: float, m: Matrix): Matrix = m + s
template `-`*(m: Matrix, s: float): Matrix = m + (-s)

proc `+=`*(m: var Matrix, s: float) =
   ## Add a matrix to a scalar in place, ``A = s+A``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = s + m[i, j]

template `-=`*(m: Matrix, s: float): Matrix = m += (-s)

proc `-`*(s: float, m: sink Matrix): Matrix =
   ## Subtract a matrix from a scalar, ``C = s-A``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = s - result[i, j]

proc `*`*(m: sink Matrix, s: float): Matrix =
   ## Multiply a matrix by a scalar, ``C = s*A``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = s * result[i, j]

proc `*=`*(m: var Matrix, s: float) =
   ## Multiply a matrix by a scalar in place, ``A = s*A``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = s * m[i, j]

template `*`*(s: float, m: Matrix): Matrix = m * s

proc `/`*(m: sink Matrix, s: float): Matrix =
   ## Divide a matrix by a scalar, ``C = A/s``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / s

template `/`*(s: float, m: Matrix): Matrix = m * (1 / s)

proc `/=`*(m: var Matrix, s: float) =
   ## Divide a matrix by a scalar in place, ``A = A/s``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = m[i, j] / s

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

proc identity*(m: int): Matrix =
   ## Generate identity matrix.
   ##
   ## ``return``: An m-by-m matrix with ones on the diagonal and zeros elsewhere.
   result = matrix(m, m)
   for i in 0 ..< m:
      result[i, i] = 1.0

template eye*(m: int): Matrix = identity(m)

proc sum*(m: Matrix): float =
   ## Sum.
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result += m[i, j]

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
   for i in 0 ..< s.len:
      let f = formatEng(s[i])
      maxLen = max(maxLen, f.len)
      result[i] = f
   for f in result.mitems:
      f = spaces(maxLen - f.len) & f

proc `$`*(m: Matrix): string =
   var formatted = newSeqOfCap[string](m.m * m.n)
   var mColj = newSeq[float](m.m)
   for j in 0 ..< m.n:
      for i in 0 ..< m.m:
         mColj[i] = m[i, j]
      formatted.add columnFormat(mColj)
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
         result.add formatted[i + j * m.m]
      if i == 0:
         result.add "⎤\n"
      elif i == m.m - 1:
         result.add "⎦"
      else:
         result.add "⎥\n"

template makeUniversal*(fname: untyped) =
   proc fname*(m: sink Matrix): Matrix =
      let len = m.m * m.n
      result = m
      for i in 0 ..< len:
         result.data[i] = fname(result.data[i])

makeUniversal(sqrt)
makeUniversal(cbrt)
makeUniversal(log10)
makeUniversal(log2)
makeUniversal(ln)
makeUniversal(exp)
makeUniversal(arccos)
makeUniversal(arcsin)
makeUniversal(arctan)
makeUniversal(cos)
makeUniversal(cosh)
makeUniversal(sin)
makeUniversal(sinh)
makeUniversal(tan)
makeUniversal(tanh)
makeUniversal(erf)
makeUniversal(erfc)
makeUniversal(lgamma)
makeUniversal(gamma)
makeUniversal(trunc)
makeUniversal(floor)
makeUniversal(ceil)
makeUniversal(degToRad)
makeUniversal(radToDeg)
