import random

type
   MatrixA* = object
      data: seq[seq[float]]
      m*, n*: int

   MatrixB* = object
      data: seq[float]
      m*, n*: int

proc `[]`*(m: MatrixA, i, j: int): float {.inline.} =
   m.data[i][j]

proc `[]`*(m: var MatrixA, i, j: int): var float {.inline.} =
   m.data[i][j]

proc `[]=`*(m: var MatrixA, i, j: int, v: float) {.inline.} =
   m.data[i][j] = v

proc rowAddr*(m: var MatrixA, i: int): ptr seq[float] =
   m.data[i].addr

proc rowUnsafeAddr*(m: MatrixA, i: int): ptr seq[float] =
   m.data[i].unsafeAddr

template newData() =
   newSeq(result.data, result.m)
   for i in 0 ..< result.m:
      newSeq(result.data[i], result.n)

proc matrixA*(m, n: int): MatrixA =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   newData()

proc randomMatrix*(m, n: Natural): MatrixA =
   const maxVal = 1000
   result.m = m
   result.n = n
   newData()
   for i in 0 ..< m:
      for j in 0 ..< n:
         result.data[i][j] = rand(maxVal).float

proc getRowPacked*(m: MatrixA): seq[float] =
   ## Make a one-dimensional row packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i * m.n + j] = m.data[i][j]

proc getColumnPacked*(m: MatrixA): seq[float] =
   ## Make a one-dimensional column packed copy of the internal array.
   newSeq(result, m.m * m.n)
   for i in 0 ..< m.m:
      for j in 0 ..< m.n:
         result[i + j * m.m] = m.data[i][j]

proc `[]`*(m: MatrixB, i, j: int): float {.inline.} =
   m.data[i * m.n + j]

proc `[]`*(m: var MatrixB, i, j: int): var float {.inline.} =
   m.data[i * m.n + j]

proc `[]=`*(m: var MatrixB, i, j: int, v: float) {.inline.} =
   m.data[i * m.n + j] = v

proc matrixB*(m, n: int): MatrixB =
   ## Construct an m-by-n matrix of zeros. 
   result.m = m
   result.n = n
   newSeq(result.data, m * n)

proc matrix*(data: seq[float], m: int): MatrixB =
   result.m = m
   result.n = if m != 0: data.len div m else: 0
   result.data = data
