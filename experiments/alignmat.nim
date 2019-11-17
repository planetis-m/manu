import random
{.passC: "-march=native -ffast-math".}
when defined(debugAsm):
   {.passC: "-fverbose-asm -masm=intel -S".}
const
   memAlign = 16

type
   Matrix* = object
      m, n: int # Row and column dimensions.
      address: pointer
      data: ptr UncheckedArray[float] # Array for internal storage of elements.

proc assumeAligned(data: ptr UncheckedArray[float]): ptr UncheckedArray[float] =
   proc builtinAssumeAligned(data: ptr UncheckedArray[float], alignment: int): ptr UncheckedArray[float] {.
      importc: "__builtin_assume_aligned", noDecl.}
   builtinAssumeAligned(data, memAlign)

{.pragma: restrict, codegenDecl: "$# __restrict $#".}
template checkBounds(cond: untyped, msg = "") =
   when compileOption("boundChecks"):
      {.line.}:
         if not cond:
            raise newException(IndexError, msg)

template allocAddress(size: int): pointer =
   alloc(size *% sizeof(float) +% (memAlign - 1))

template alignData(address: pointer): ptr UncheckedArray[float] =
   let aligned {.restrict.} = block:
      echo cast[int](address)
      let remainder = cast[ByteAddress](address) and (memAlign - 1)
      echo remainder
      if remainder == 0:
         assumeAligned cast[ptr UncheckedArray[float]](address)
      else:
         let offset = memAlign -% remainder
         echo offset
         echo cast[ByteAddress](address) +% offset
         assumeAligned cast[ptr UncheckedArray[float]](cast[ByteAddress](address) +% offset)
   aligned

proc `=destroy`*(m: var Matrix) =
   if m.address != nil:
      dealloc(m.address)
      m.address = nil
      m.data = nil
      m.m = 0
      m.n = 0

proc `=sink`*(a: var Matrix; b: Matrix) =
   `=destroy`(a)
   a.address = b.address
   a.data = b.data
   a.m = b.m
   a.n = b.n

proc `=`*(a: var Matrix; b: Matrix) =
   if a.address != b.address:
      `=destroy`(a)
      a.m = b.m
      a.n = b.n
      if b.address != nil:
         let len = b.m * b.n
         a.address = allocAddress(len)
         a.data = alignData(a.address)
         copyMem(a.address, b.address, len * sizeof(float) + memAlign - 1)

proc matrix*(m, n: int): Matrix =
   ## Construct an m-by-n matrix of zeros.
   result.m = m
   result.n = n
   let len = m * n
   result.address = allocAddress(len)
   result.data = alignData(result.address)
   for i in 0 ..< len:
      result.data[i] = 0.0

proc matrix*(m, n: int, s: float): Matrix =
   ## Construct an m-by-n constant matrix.
   result.m = m
   result.n = n
   let len = m * n
   result.address = allocAddress(len)
   result.data = alignData(result.address)
   for i in 0 ..< len:
      result.data[i] = s

proc matrix*(data: seq[seq[float]]): Matrix =
   ## Construct a matrix from a 2-D array.
   result.m = data.len
   result.n = data[0].len
   for i in 0 ..< result.m:
      assert(data[i].len == result.n, "All rows must have the same length.")
   result.address = allocAddress(result.m * result.n)
   result.data = alignData(result.address)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result.data[i * result.n + j] = data[i][j]

proc matrix*(data: seq[seq[float]], m, n: int): Matrix =
   ## Construct a matrix quickly without checking arguments.
   result.m = m
   result.n = n
   let len = m * n
   result.address = allocAddress(len)
   result.data = alignData(result.address)
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
   result.address = allocAddress(data.len)
   result.data = alignData(result.address)
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i * n + j] = data[i + j * m]

proc randMatrix*(m, n: int): Matrix =
   ## Generate matrix with random elements.
   ##
   ## ``return``: an m-by-n matrix with uniformly distributed random elements.
   result.m = m
   result.n = n
   let len = m * n
   result.address = allocAddress(len)
   result.data = alignData(result.address)
   for i in 0 ..< len:
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
   let data {.restrict.} = assumeAligned m.data
   data[i * m.n + j]

proc `[]`*(m: var Matrix, i, j: int): var float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   let data {.restrict.} = assumeAligned m.data
   data[i * m.n + j]

proc `[]=`*(m: var Matrix, i, j: int, s: float) {.inline.} =
   ## Set a single element.
   checkBounds(i >= 0 and i < m.m)
   checkBounds(j >= 0 and j < m.n)
   let data {.restrict.} = assumeAligned m.data
   data[i * m.n + j] = s

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
   let b = a
   let c = a *. b
   echo c[3, 4]


main()
