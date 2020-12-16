import random, strutils
{.passC: "-march=native -ffast-math".}
when defined(debugAsm):
   {.passC: "-fverbose-asm -masm=intel -S".}

type
   Matrix* = object
      m*, n*: int # Row and column dimensions.
      data: ptr UncheckedArray[float] # Array for internal storage of elements.
type All* = object
template checkBounds(cond: untyped, msg = "") =
   when compileOption("boundChecks"):
      {.line.}:
         if not cond:
            raise newException(IndexDefect, msg)

template createData(size): ptr UncheckedArray[float] =
   cast[ptr UncheckedArray[float]](alloc0(size * sizeof(float)))

proc `=destroy`*(m: var Matrix) =
   if m.data != nil:
      dealloc(m.data)

proc `=copy`*(a: var Matrix; b: Matrix) =
   if a.data != b.data:
      `=destroy`(a)
      wasMoved(a)
      a.m = b.m
      a.n = b.n
      if b.data != nil:
         let len = b.m * b.n
         a.data = createData(len)
         copyMem(a.data, b.data, len * sizeof(float))

template printData(m) =
   echo cast[ByteAddress](m.data).toHex

proc matrix*(m, n: int): Matrix =
   ## Construct an m-by-n matrix of zeros.
   result.m = m
   result.n = n
   result.data = createData(m * n)

proc matrix*(m, n: int, s: float): Matrix =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   let len = m * n
   result.data = createData(len)
   for i in 0 ..< len:
      result.data[i] = s

proc matrix*(n: int, data: seq[float]): Matrix =
   ## Construct a matrix from a one-dimensional packed array.
   ##
   ## parameter ``data``: one-dimensional array of SomeFloat, packed by rows.
   ## Array length must be a multiple of ``n``.
   let m = if n != 0: data.len div n else: 0
   assert(m * n == data.len, "Array length must be a multiple of n.")
   result.m = m
   result.n = n
   result.data = createData(data.len)
   #for i in 0 ..< data.len:
      #result.data[i] = data[i]
   copyMem(result.data, data[0].unsafeAddr, data.len * sizeof(float))

proc randMatrix*(m, n: int): Matrix =
   ## Generate matrix with random elements.
   ##
   ## ``return``: an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   let len = m * n
   result.data = createData(len)
   for i in 0 ..< len:
      result.data[i] = rand(1.0)

proc getColumnPacked*(m: Matrix): seq[float] =
   ## Make a one-dimensional column packed copy of the internal array.
   result = newSeq[float](m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i + j * m.m] = m.data[i * m.n + j]

proc getRowPacked*(m: Matrix): seq[float] =
   ## Copy the internal one-dimensional row packed array.
   result = newSeq[float](m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i * m.n + j] = m.data[i * m.n + j]

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

proc `[]`*[U, V](m: Matrix, r: HSlice[U, V], c: typedesc[All]): lent Matrix =
   ## Get a submatrix, all columns
   ## ``m[i0 .. i1, 0 .. ^1]``
   let ra = m.m ^^ r.a
   let rb = m.m ^^ r.b
   checkBounds(ra >= 0 and rb < m.m, "Submatrix dimensions")
   printData(m)
   result = m
   printData(result)
   result.m = rb - ra + 1
   echo cast[ByteAddress](result.data[ra].addr).toHex
   result.data = cast[ptr UncheckedArray[float]](result.data[ra].addr)
   printData(result)

proc `+`*(a: sink Matrix; b: Matrix): Matrix =
   ## ``C = A + B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[i, j]

proc `-`*(m: sink Matrix): Matrix =
   ## Unary minus
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = -result[i, j]

proc `-`*(a: sink Matrix; b: Matrix): Matrix =
   ## ``C = A - B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[i, j]

proc `*.`*(a: sink Matrix; b: Matrix): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[i, j]

proc main =
   let a = matrix(5, 5, 4.0)
   let b = matrix(2, @[0.0, 0, 0, 1, 1, 0, 1, 1])
   let c = b[1..1, All]
   printData(c)


main()
