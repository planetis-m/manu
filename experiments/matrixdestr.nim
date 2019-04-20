type
   Matrix* = object
      m*, n*: int # Row and column dimensions.
      data: ptr UncheckedArray[float] # Array for internal storage of elements.

proc `=destroy`*(m: var Matrix) =
   if m.data != nil:
      dealloc(m.data)
      m.data = nil
      m.m = 0
      m.n = 0

proc `=sink`*(a: var Matrix; b: Matrix) =
   if a.data != nil and a.data != b.data:
      dealloc(a.data)
   a.data = b.data
   a.m = b.m
   a.n = b.n

proc `=`*(a: var Matrix; b: Matrix) =
   if a.data != nil and a.data != b.data:
      dealloc(a.data)
      a.data = nil
   a.m = b.m
   a.n = b.n
   if b.data != nil:
      a.data = cast[type(a.data)](alloc(a.m * a.n))
      copyMem(a.data, b.data, b.m * b.n)

proc matrix*(m, n: int): Matrix =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   result.data = cast[type(result.data)](alloc(m * n))

proc matrix*(m, n: int, s: float): Matrix =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   result.data = cast[type(result.data)](alloc(m * n))
   for i in 0 ..< m * n:
      result.data[i] = s

proc matrix*(data: seq[seq[float]]): Matrix =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
   for i in 0 ..< result.m:
      assert(data[i].len == result.n, "All rows must have the same length.")
   result.data = cast[type(result.data)](alloc(result.m * result.n))
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result.data[i * result.n + j] = data[i][j]

proc matrix*(data: seq[seq[float]], m, n: int): Matrix =
   ## Construct a matrix quickly without checking arguments.
   result.m = m
   result.n = n
   result.data = cast[type(result.data)](alloc(m * n))
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
   result.data = cast[type(result.data)](alloc(data.len))
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i + j * m]

proc randMatrix*(m, n: int): Matrix =
   ## Generate matrix with random elements.
   ##
   ## ``return``: an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   result.data = cast[type(result.data)](alloc(m * n))
   for i in 0 ..< m * n:
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
   result = newSeq[float](m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i * m.n + j] = m.data[i * m.n + j]
