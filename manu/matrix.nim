import math, random, strutils

template checkBounds(cond: untyped, msg = "") =
   when compileOption("boundChecks"):
      {.line.}:
         if not cond:
            raise newException(IndexDefect, msg)

type
   Matrix*[T: SomeFloat] = object
      m, n: int # Row and column dimensions.
      data: ptr UncheckedArray[T] # Array for internal storage of elements.

template createData[T](size): ptr UncheckedArray[T] =
   cast[ptr UncheckedArray[T]](allocShared(size * sizeof(T)))

template createData0[T](size): ptr UncheckedArray[T] =
   cast[ptr UncheckedArray[T]](allocShared0(size * sizeof(T)))

proc `=destroy`*[T](m: var Matrix[T]) =
   if m.data != nil:
      dealloc(m.data)

proc `=`*[T](a: var Matrix[T]; b: Matrix[T]) =
   if a.data != b.data:
      `=destroy`(a)
      wasMoved(a)
      a.m = b.m
      a.n = b.n
      if b.data != nil:
         let len = b.m * b.n
         a.data = createData[T](len)
         copyMem(a.data, b.data, len * sizeof(T))

type
   Matrix64* = Matrix[float64]
      ## Alias for a ``Matrix`` of 64-bit floats.
   Matrix32* = Matrix[float32]
      ## Alias for a ``Matrix`` of 32-bit floats.
   ColVector*[T: SomeFloat] = distinct Matrix[T]
   RowVector*[T: SomeFloat] = distinct Matrix[T]
   ColVector64* = ColVector[float64]
   RowVector64* = RowVector[float64]
   ColVector32* = ColVector[float32]
   RowVector32* = RowVector[float32]

proc matrix*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
   ## Construct an m-by-n matrix of zeros.
   result.m = m
   result.n = n
   result.data = createData0[T](m * n)

proc matrixUninit*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
   ## Construct an m-by-n matrix. Note that the matrix will be uninitialized.
   result.m = m
   result.n = n
   result.data = createData[T](m * n)

proc matrix*[T: SomeFloat](m, n: int, s: T): Matrix[T] =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   let len = m * n
   result.data = createData[T](len)
   for i in 0 ..< len:
      result.data[i] = s

template ones*[T](m, n: int): Matrix[T] = matrix[T](m, n, T(1.0))
template zeros*[T](m, n: int): Matrix[T] = matrix[T](m, n)
template ones32*(m, n: int): Matrix32 = matrix(m, n, 1.0'f32)
template zeros32*(m, n: int): Matrix32 = matrix[float32](m, n)
template ones64*(m, n: int): Matrix64 = matrix(m, n, 1.0)
template zeros64*(m, n: int): Matrix64 = matrix[float64](m, n)

proc matrix*[T: SomeFloat](data: seq[seq[T]]): Matrix[T] =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
   for i in 0 ..< result.m:
      assert(data[i].len == result.n, "All rows must have the same length.")
   result.data = createData[T](result.m * result.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result.data[i * result.n + j] = data[i][j]

proc matrix*[T: SomeFloat](data: seq[seq[T]], m, n: int): Matrix[T] =
   ## Construct a matrix quickly without checking arguments.
   result.m = m
   result.n = n
   result.data = createData[T](m * n)
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i][j]

proc matrix*[T: SomeFloat](data: seq[T], m: int): Matrix[T] =
   ## Construct a matrix from a one-dimensional packed array.
   ##
   ## parameter ``data``: one-dimensional array of float, packed by columns (ala Fortran).
   ## Array length must be a multiple of ``m``.
   let n = if m != 0: data.len div m else: 0
   assert(m * n == data.len, "Array length must be a multiple of m.")
   result.m = m
   result.n = n
   result.data = createData[T](data.len)
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i + j * m]

proc matrix*[T: SomeFloat](n: int, data: seq[T]): Matrix[T] =
   ## Construct a matrix from a one-dimensional packed array.
   ##
   ## parameter ``data``: one-dimensional array of SomeFloat, packed by rows.
   ## Array length must be a multiple of ``n``.
   let m = if n != 0: data.len div n else: 0
   assert(m * n == data.len, "Array length must be a multiple of n.")
   result.m = m
   result.n = n
   result.data = createData[T](data.len)
   for i in 0 ..< data.len:
      result.data[i] = data[i]
#    copyMem(result.data, data[0].unsafeAddr, data.len * sizeof(T))

proc randMatrix*[T: SomeFloat](m, n: int, max: T): Matrix[T] =
   ## Generate matrix with random elements.
   ##
   ## ``return``: an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   let len = m * n
   result.data = createData[T](len)
   for i in 0 ..< len:
      result.data[i] = rand(max)

proc randMatrix*[T: SomeFloat](m, n: int, x: Slice[T]): Matrix[T] =
   ## Generate matrix with random elements.
   ##
   ## ``return``: an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   let len = m * n
   result.data = createData[T](len)
   for i in 0 ..< len:
      result.data[i] = rand(x)

template randMatrix*[T](m, n: int): Matrix[T] = randMatrix[T](m, n, T(1.0))
template randMatrix32*(m, n: int): Matrix32 = randMatrix(m, n, 1.0'f32)
template randMatrix64*(m, n: int): Matrix64 = randMatrix(m, n, 1.0)

when not declared(gauss):
   proc gauss*(r: var Rand; mu = 0.0; sigma = 1.0): float =
      # Ratio of uniforms method for normal
      # http://www2.econ.osaka-u.ac.jp/~tanizaki/class/2013/econome3/13.pdf
      const K = sqrt(2 / E)
      var
         a = 0.0
         b = 0.0
      while true:
         a = rand(r, 1.0)
         b = (2.0 * rand(r, 1.0) - 1.0) * K
         if  b * b <= -4.0 * a * a * ln(a): break
      result = mu + sigma * (b / a)

proc randNMatrix*[T: SomeFloat](m, n: int, mu, sigma: T): Matrix[T] =
   ## Normal distribution
   ##
   ## parameter ``mu``: the mean
   ## parameter ``sigma``: the standard deviation
   ## ``return``: an m-by-n matrix with normally distributed random elements.
   result.m = m
   result.n = n
   let len = m * n
   result.data = createData[T](len)
   for i in 0 ..< len:
      result.data[i] = gauss(mu, sigma).T

template randNMatrix*[T](m, n: int): Matrix[T] = randNMatrix[T](m, n, T(0.0), T(1.0))
template randNMatrix32*(m, n: int): Matrix32 = randNMatrix(m, n, 0.0'f32, 1.0'f32)
template randNMatrix64*(m, n: int): Matrix64 = randNMatrix(m, n, 0.0, 1.0)

proc getArray*[T](m: Matrix[T]): seq[seq[T]] =
   ## Make a two-dimensional array copy of the internal array.
   result = newSeq[seq[T]](m.m)
   for i in 0 ..< m.m:
      result[i] = newSeq[T](m.n)
      for j in 0 ..< m.n:
         result[i][j] = m.data[i * m.n + j]

proc getColumnPacked*[T](m: Matrix[T]): seq[T] =
   ## Make a one-dimensional column packed copy of the internal array.
   result = newSeq[T](m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i + j * m.m] = m.data[i * m.n + j]

proc getRowPacked*[T](m: Matrix[T]): seq[T] {.inline.} =
   ## Copy the internal one-dimensional row packed array.
   result = newSeq[T](m.m * m.n)
   for i in 0 ..< result.len:
      result[i] = m.data[i]

proc dim*[T](m: Matrix[T]): (int, int) {.inline.} =
   ## Get (row, column) dimensions tuple.
   (m.m, m.n)

proc m*[T](m: Matrix[T]): int {.inline.} =
   ## Get row dimension.
   m.m

proc n*[T](m: Matrix[T]): int {.inline.} =
   ## Get column dimension.
   m.n

proc rowDimension*[T](m: Matrix[T]): int {.inline.} = m.m
proc columnDimension*[T](m: Matrix[T]): int {.inline.} = m.n

proc `[]`*[T](m: Matrix[T], i, j: int): lent T {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   result = m.data[i * m.n + j]

proc `[]`*[T](m: var Matrix[T], i, j: int): var T {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   m.data[i * m.n + j]

proc `[]=`*[T](m: var Matrix[T], i, j: int, s: T) {.inline.} =
   ## Set a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   m.data[i * m.n + j] = s

proc `[]`*[T](v: ColVector[T], i: int): lent T {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < Matrix[T](v).m)
   Matrix[T](v).data[i]

proc `[]`*[T](v: var ColVector[T], i: int): var T {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < Matrix[T](v).m)
   Matrix[T](v).data[i]

proc `[]=`*[T](v: var ColVector[T], i: int, s: T) {.inline.} =
   ## Set a single element.
   checkBounds(i >= 0 and i < Matrix[T](v).m)
   Matrix[T](v).data[i] = s

proc `[]`*[T](v: RowVector[T], j: int): lent T {.inline.} =
   ## Get a single element.
   checkBounds(j >= 0 and j < Matrix[T](v).n)
   Matrix[T](v).data[j]

proc `[]`*[T](v: var RowVector[T], j: int): var T {.inline.} =
   ## Get a single element.
   checkBounds(j >= 0 and j < Matrix[T](v).n)
   Matrix[T](v).data[j]

proc `[]=`*[T](v: var RowVector[T], j: int, s: T) {.inline.} =
   ## Set a single element.
   checkBounds(j >= 0 and j < Matrix[T](v).n)
   Matrix[T](v).data[j] = s

template `^^`(dim, i: untyped): untyped =
  (when i is BackwardsIndex: dim - int(i) else: int(i))

proc `[]`*[T, U, V, W, X](m: Matrix[T], r: HSlice[U, V], c: HSlice[W, X]): Matrix[T] =
   ## Get a submatrix,
   ## ``m[i0 .. i1, j0 .. j1]``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(ra >= 0 and rb < m.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(ca >= 0 and cb < m.n, "Submatrix dimensions")
   result = matrixUninit[T](rb - ra + 1, cb - ca + 1)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = m[i + ra, j + ca]

proc `[]`*[T](m: Matrix[T], r, c: openarray[SomeInteger]): Matrix[T] =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result = matrixUninit[T](r.len, c.len)
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         result[i, j] = m[r[i], c[j]]

proc `[]`*[T, U, V](m: Matrix[T], r: HSlice[U, V], c: openarray[SomeInteger]): Matrix[T] =
   ## Get a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]]``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(ra >= 0 and rb < m.m, "Submatrix dimensions")
   checkBounds(c.len <= m.n, "Submatrix dimensions")
   result = matrixUninit[T](rb - ra + 1, c.len)
   for i in 0 ..< result.m:
      for j in 0 ..< c.len:
         result[i, j] = m[i + ra, c[j]]

proc `[]`*[T, U, V](m: Matrix[T], r: openarray[SomeInteger], c: HSlice[U, V]): Matrix[T] =
   ## Get a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1]``
   checkBounds(r.len <= m.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(ca >= 0 and cb < m.n, "Submatrix dimensions")
   result = matrixUninit[T](r.len, cb - ca + 1)
   for i in 0 ..< r.len:
      for j in 0 ..< result.n:
         result[i, j] = m[r[i], j + ca]

proc `[]=`*[T, U, V, W, X](m: var Matrix[T], r: HSlice[U, V], c: HSlice[W, X], a: Matrix[T]) =
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

proc `[]=`*[T](m: var Matrix[T], r, c: openarray[SomeInteger], a: Matrix[T]) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], [1, 2, 3, 4]] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in 0 ..< c.len:
         m[r[i], c[j]] = a[i, j]

proc `[]=`*[T, U, V](m: var Matrix[T], r: HSlice[U, V], c: openarray[SomeInteger], a: Matrix[T]) =
   ## Set a submatrix,
   ## ``m[i0 .. i1, [0, 2, 3, 4]] = a``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(rb - ra + 1 == a.m, "Submatrix dimensions")
   checkBounds(c.len == a.n, "Submatrix dimensions")
   for i in 0 ..< a.m:
      for j in 0 ..< c.len:
         m[i + ra, c[j]] = a[i, j]

proc `[]=`*[T, U, V](m: var Matrix[T], r: openarray[SomeInteger], c: HSlice[U, V], a: Matrix[T]) =
   ## Set a submatrix,
   ## ``m[[0, 2, 3, 4], j0 .. j1] = a``
   checkBounds(r.len == a.m, "Submatrix dimensions")
   let ca = m.n ^^ c.a
   let cb = m.n ^^ c.b
   checkBounds(cb - ca + 1 == a.n, "Submatrix dimensions")
   for i in 0 ..< r.len:
      for j in 0 ..< a.n:
         m[r[i], j + ca] = a[i, j]

proc `-`*[T](m: sink Matrix[T]): Matrix[T] =
   ## Unary minus
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = -result[i, j]

proc `+`*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
   ## ``C = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[i, j]

proc `+=`*[T](a: var Matrix[T], b: Matrix[T]) =
   ## ``A = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[i, j]

proc `-`*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
   ## ``C = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[i, j]

proc `-=`*[T](a: var Matrix[T], b: Matrix[T]) =
   ## ``A = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[i, j]

proc `*.`*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[i, j]

proc `*.=`*[T](a: var Matrix[T], b: Matrix[T]) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[i, j]

proc `/.`*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
   ## Element-by-element right division, ``C = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / b[i, j]

proc `/.=`*[T](a: var Matrix[T], b: Matrix[T]) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[i, j]

proc `\.`*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[i, j] / result[i, j]

proc `\.=`*[T](a: var Matrix[T], b: Matrix[T]) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[i, j] / a[i, j]

proc `+`*[T](m: sink Matrix[T], s: T): Matrix[T] =
   ## Add a matrix to a scalar, ``C = s+A``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = s + result[i, j]

template `+`*[T](s: T, m: Matrix[T]): Matrix[T] = m + s
template `-`*[T](m: Matrix[T], s: T): Matrix[T] = m + (-s)

proc `+=`*[T](m: var Matrix[T], s: T) =
   ## Add a matrix to a scalar in place, ``A = s+A``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = s + m[i, j]

template `-=`*[T](m: Matrix[T], s: T): Matrix[T] = m += (-s)

proc `-`*[T](s: T, m: sink Matrix[T]): Matrix[T] =
   ## Subtract a matrix from a scalar, ``C = s-A``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = s - result[i, j]

proc `*`*[T](m: sink Matrix[T], s: T): Matrix[T] =
   ## Multiply a matrix by a scalar, ``C = s*A``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = s * result[i, j]

proc `*=`*[T](m: var Matrix[T], s: T) =
   ## Multiply a matrix by a scalar in place, ``A = s*A``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = s * m[i, j]

template `*`*[T](s: T, m: Matrix[T]): Matrix[T] = m * s

proc `/`*[T](m: sink Matrix[T], s: T): Matrix[T] =
   ## Divide a matrix by a scalar, ``C = A/s``
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / s

template `/`*[T](s: T, m: Matrix[T]): Matrix[T] = m * (1 / s)

proc `/=`*[T](m: var Matrix[T], s: T) =
   ## Divide a matrix by a scalar in place, ``A = A/s``
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         m[i, j] = m[i, j] / s

proc `*`*[T](a, b: Matrix[T]): Matrix[T] =
   ## Linear algebraic matrix multiplication, ``A * B``
   assert(b.m == a.n, "Matrix inner dimensions must agree.")
   result = matrixUninit[T](a.m, b.n)
   var bColj = newSeq[T](b.m)
   for j in 0 ..< b.n:
      for k in 0 ..< a.n:
         bColj[k] = b[k, j]
      for i in 0 ..< a.m:
         var s = 0.0
         for k in 0 ..< a.n:
            s += a[i, k] * bColj[k]
         result[i, j] = s

proc `+`*[T](a: sink Matrix[T], b: ColVector[T]): Matrix[T] =
   ## ``C = A + B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[i]

proc `+`*[T](a: sink Matrix[T], b: RowVector[T]): Matrix[T] =
   ## ``C = A + B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[j]

proc `+`*[T](a: ColVector[T], b: RowVector[T]): Matrix[T] =
   ## ``C = A + B``, ``a`` and ``b`` are broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](a).n == 1, "Matrices are not broadcastable.")
   result = matrix(Matrix[T](a).m, Matrix[T](b).n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] + b[j]

template `+`*[T](b: RowVector[T], a: ColVector[T]): Matrix[T] = a + b
template `+`*[T](b: ColVector[T], a: Matrix[T]): Matrix[T] = a + b
template `+`*[T](b: RowVector[T], a: Matrix[T]): Matrix[T] = a + b

proc `+=`*[T](a: var Matrix[T], b: ColVector[T]) =
   ## ``A = A + B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[i]

proc `+=`*[T](a: var Matrix[T], b: RowVector[T]) =
   ## ``A = A + B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[j]

proc `-`*[T](a: sink Matrix[T], b: ColVector[T]): Matrix[T] =
   ## ``C = A - B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[i]

proc `-`*[T](a: sink Matrix[T], b: RowVector[T]): Matrix[T] =
   ## ``C = A - B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[j]

proc `-`*[T](a: ColVector[T], b: RowVector[T]): Matrix[T] =
   ## ``C = A - B``, ``a`` and ``b`` are broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](a).n == 1, "Matrices are not broadcastable.")
   result = matrix(Matrix[T](a).m, Matrix[T](b).n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] - b[j]

proc `-`*[T](a: ColVector[T], b: sink Matrix[T]): Matrix[T] =
   ## ``C = A - B``, ``b`` is broadcasted
   assert(Matrix[T](a).m == b.m and Matrix[T](a).n == 1, "Matrix-vector dimensions must agree.")
   result = b
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] - result[i, j]

proc `-`*[T](a: RowVector[T], b: sink Matrix[T]): Matrix[T] =
   ## ``C = A - B``, ``b`` is broadcasted
   assert(Matrix[T](a).m == 1 and Matrix[T](a).n == b.n, "Matrix-vector dimensions must agree.")
   result = b
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[j] - result[i, j]

proc `-`*[T](a: RowVector[T], b: ColVector[T]): Matrix[T] =
   ## ``C = A - B``, ``a`` and ``b`` are broadcasted
   assert(Matrix[T](a).m == 1 and Matrix[T](b).n == 1, "Matrices are not broadcastable.")
   result = matrix(Matrix[T](b).m, Matrix[T](a).n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[j] - b[i]

proc `-=`*[T](a: var Matrix[T], b: ColVector[T]) =
   ## ``A = A - B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[i]

proc `-=`*[T](a: var Matrix[T], b: RowVector[T]) =
   ## ``A = A - B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[j]

proc `*.`*[T](a: sink Matrix[T], b: ColVector[T]): Matrix[T] =
   ## Element-by-element multiplication, ``C = A.*B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[i]

proc `*.`*[T](a: sink Matrix[T], b: RowVector[T]): Matrix[T] =
   ## Element-by-element multiplication, ``C = A.*B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[j]

proc `*.`*[T](a: ColVector[T], b: RowVector[T]): Matrix[T] =
   ## Element-by-element multiplication, ``C = A.*B``, ``a`` and ``b`` are broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](a).n == 1, "Matrices are not broadcastable.")
   result = matrix(Matrix[T](a).m, Matrix[T](b).n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] * b[j]

template `*.`*[T](b: RowVector[T], a: ColVector[T]): Matrix[T] = a *. b
template `*.`*[T](b: ColVector[T], a: Matrix[T]): Matrix[T] = a *. b
template `*.`*[T](b: RowVector[T], a: Matrix[T]): Matrix[T] = a *. b

proc `*.=`*[T](a: var Matrix[T], b: ColVector[T]) =
   ## Element-by-element multiplication in place, ``A = A.*B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[i]

proc `*.=`*[T](a: var Matrix[T], b: RowVector[T]) =
   ## Element-by-element multiplication in place, ``A = A.*B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[j]

proc `/.`*[T](a: sink Matrix[T], b: ColVector[T]): Matrix[T] =
   ## Element-by-element right division, ``C = A./B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / b[i]

proc `/.`*[T](a: sink Matrix[T], b: RowVector[T]): Matrix[T] =
   ## Element-by-element right division, ``C = A./B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / b[j]

proc `/.`*[T](a: ColVector[T], b: RowVector[T]): Matrix[T] =
   ## Element-by-element right division, ``C = A./B``, ``a`` and ``b`` are broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](a).n == 1, "Matrices are not broadcastable.")
   result = matrix(Matrix[T](a).m, Matrix[T](b).n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] / b[j]

proc `/.=`*[T](a: var Matrix[T], b: ColVector[T]) =
   ## Element-by-element right division in place, ``A = A./B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[i]

proc `/.=`*[T](a: var Matrix[T], b: RowVector[T]) =
   ## Element-by-element right division in place, ``A = A./B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[j]

proc `\.`*[T](a: sink Matrix[T], b: ColVector[T]): Matrix[T] =
   ## Element-by-element left division, ``C = A.\B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[i] / result[i, j]

proc `\.`*[T](a: sink Matrix[T], b: RowVector[T]): Matrix[T] =
   ## Element-by-element left division, ``C = A.\B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[j] / result[i, j]

proc `\.`*[T](a: ColVector[T], b: RowVector[T]): Matrix[T] =
   ## Element-by-element left division, ``C = A.\B``, ``a`` and ``b`` are broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](a).n == 1, "Matrices are not broadcastable.")
   result = matrix(Matrix[T](a).m, Matrix[T](b).n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[j] / a[i]

proc `\.=`*[T](a: var Matrix[T], b: ColVector[T]) =
   ## Element-by-element left division in place, ``A = A.\B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == a.m and Matrix[T](b).n == 1, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[i] / a[i, j]

proc `\.=`*[T](a: var Matrix[T], b: RowVector[T]) =
   ## Element-by-element left division in place, ``A = A.\B``, ``b`` is broadcasted
   assert(Matrix[T](b).m == 1 and Matrix[T](b).n == a.n, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[j] / a[i, j]

template `/.`*[T](b: RowVector[T], a: ColVector[T]): Matrix[T] = a \. b
template `/.`*[T](b: ColVector[T], a: Matrix[T]): Matrix[T] = a \. b
template `/.`*[T](b: RowVector[T], a: Matrix[T]): Matrix[T] = a \. b

template `\.`*[T](b: RowVector[T], a: ColVector[T]): Matrix[T] = a /. b
template `\.`*[T](b: ColVector[T], a: Matrix[T]): Matrix[T] = a /. b
template `\.`*[T](b: RowVector[T], a: Matrix[T]): Matrix[T] = a /. b

proc transpose*[T](m: Matrix[T]): Matrix[T] =
   ## Matrix transpose
   result = matrixUninit[T](m.n, m.m)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[j, i] = m[i, j]

proc identity*[T: SomeFloat](m: int): Matrix[T] =
   ## Generate identity matrix.
   ##
   ## ``return``: An m-by-m matrix with ones on the diagonal and zeros elsewhere.
   result = matrix[T](m, m)
   for i in 0 ..< m:
      result[i, i] = T(1.0)

template eye*[T](m: int): Matrix[T] = identity[T](m)
template eye32*(m: int): Matrix32 = identity[float32](m)
template eye64*(m: int): Matrix64 = identity[float64](m)

proc sum*[T](m: Matrix[T]): T =
   ## Sum of all elements.
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result += m[i, j]

proc sumColumns*[T](m: Matrix[T]): Matrix[T] =
   ## Sum of matrix columns.
   ##
   ## ``return``: An 1-by-n matrix with sum of the elements in each column.
   result = matrixUninit[T](1, m.n)
   for j in 0 ..< m.n:
      var s = T(0.0)
      for i in 0 ..< m.m:
         s += m[i, j]
      RowVector[T](result)[j] = s

proc sumRows*[T](m: Matrix[T]): Matrix[T] =
   ## Sum of matrix rows.
   ##
   ## ``return``: An m-by-1 matrix with sum of the elements in each row.
   result = matrixUninit[T](m.m, 1)
   for i in 0 ..< m.m:
      var s = T(0.0)
      for j in 0 ..< m.n:
         s += m[i, j]
      ColVector[T](result)[i] = s

proc norm1*[T](m: Matrix[T]): T =
   ## One norm.
   ##
   ## ``return``: maximum column sum
   for j in 0 ..< m.n:
      var s = T(0.0)
      for i in 0 ..< m.m:
         s += abs(m[i, j])
      result = max(result, s)

proc normInf*[T](m: Matrix[T]): T =
   ## Infinity norm.
   ##
   ## ``return``: maximum row sum
   for i in 0 ..< m.m:
      var s = T(0.0)
      for j in 0 ..< m.n:
         s += abs(m[i, j])
      result = max(result, s)

proc normF*[T](m: Matrix[T]): T =
   ## Frobenius norm.
   ##
   ## ``return``: sqrt of sum of squares of all elements.
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result = hypot(result, m[i, j])

proc trace*[T](m: Matrix[T]): T =
   ## Matrix trace.
   ##
   ## ``return``: the sum of the diagonal elements
   for i in 0 ..< min(m.m, m.n):
      result += m[i, i]

proc columnFormat[T](s: seq[T]): seq[string] =
   result = newSeq[string](s.len)
   var maxLen = 0
   for i in 0 ..< s.len:
      let f = formatEng(s[i])
      maxLen = max(maxLen, f.len)
      result[i] = f
   for f in result.mitems:
      f = spaces(maxLen - f.len) & f

proc `$`*[T](m: Matrix[T]): string =
   var formatted = newSeqOfCap[string](m.m * m.n)
   var mColj = newSeq[T](m.m)
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
   proc fname*[T](m: sink Matrix[T]): Matrix[T] =
      result = m
      for i in 0 ..< result.m:
         for j in 0 ..< result.n:
            result[i, j] = fname(result[i, j])

template makeUniversalBinary*(fname: untyped) =
   proc fname*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
      assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
      result = a
      for i in 0 ..< result.m:
         for j in 0 ..< result.n:
            result[i, j] = fname(result[i, j], b[i, j])

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
