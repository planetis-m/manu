# const
#   memAlign {.intdefine.} = 32

type
  # MatrixStorageOrder* = enum
  #   RowMajorStorage, ColMajorStorage

  MatPayloadBase = object
    counter: int

  MatPayload = object
    counter: int
    data: UncheckedArray[float]

  Matrix* = object
    len, offset, stride: int
    p: ptr MatPayload # can be nil if len == 0.

template checkBounds(cond: untyped, msg = "") =
  when compileOption("boundChecks"):
    {.line.}:
      if not cond:
        raise newException(IndexDefect, msg)

template contentSize(cap): int = cap*sizeof(float) + sizeof(MatPayloadBase)

template frees(s) =
  when compileOption("threads"):
    deallocShared(s.p)
  else:
    dealloc(s.p)

proc `=destroy`*(x: var Matrix) =
  if x.p != nil:
    if x.p.counter == 0:
      frees(x)
    else:
      dec x.p.counter

template dups(a, b) =
  if b.p != nil:
    inc b.p.counter
  a.p = b.p
  a.len = b.len

proc `=dup`*(b: Matrix): Matrix =
  dups(result, b)

proc `=copy`*(a: var Matrix, b: Matrix) =
  `=destroy`(a)
  dups(a, b)

template allocs(s) =
  when compileOption("threads"):
    let p {.inject.} = cast[ptr MatPayload](allocShared(s))
  else:
    let p {.inject.} = cast[ptr MatPayload](alloc(s))

template allocs0(s) =
  when compileOption("threads"):
    let p {.inject.} = cast[ptr MatPayload](allocShared0(s))
  else:
    let p {.inject.} = cast[ptr MatPayload](alloc0(s))

proc deepCopy*(y: Matrix): Matrix =
  if y.len <= 0:
    result = Matrix(len: 0, offset: 0, stride: 0, p: nil)
  else:
    allocs(contentSize(y.len))
    p.counter = 0
    copyMem(addr p.data[0], addr y.p.data[y.offset], y.len * sizeof(T))
    result = Matrix(len: y.len, offset: 0, stride: y.stride, p: p)

proc matrix*(m, n: int): Matrix {.inline.} =
  ## Construct an m-by-n matrix of zeros.
  let len = m * n
  allocs0(contentSize(len))
  result = Matrix(len: len, offset: 0, stride: n, p: p)

proc matrix*(data: seq[float], m: int): Matrix =
  ## Construct a matrix from a one-dimensional packed array.
  ##
  ## parameter ``data``: one-dimensional array of float, packed by columns (ala Fortran).
  ## Array length must be a multiple of ``m``.
  let n = if m != 0: data.len div m else: 0
  let len = m * n
  assert(len == data.len, "Array length must be a multiple of m.")
  allocs(contentSize(len))
  # for i in 0 ..< m:
  #   for j in 0 ..< n:
  #     p.data[i * n + j] = data[i + j * m]
  var rx = 0
  for bx in countup(0, m-1):
    for ax in countup(0, high(data), m):
      p.data[rx] = data[ax + bx]
      inc rx
  result = Matrix(len: len, offset: 0, stride: n, p: p)

proc m*(m: Matrix): int {.inline.} =
  ## Get row dimension.
  m.len div m.stride

proc n*(m: Matrix): int {.inline.} =
  ## Get column dimension.
  m.stride

proc `[]`*(m: Matrix, i, j: int): float {.inline.} =
  ## Get a single element.
  checkBounds(i >= 0 and i < m.len div m.stride) # lol
  checkBounds(j >= 0 and j < m.stride)
  result = m.p.data[m.offset + i * m.stride + j]

import std/strutils

proc `$`(m: Matrix): string =
  result = ""
  for e in countup(0, m.len-1, m.stride):
    if e > 0: result.add "\n"
    var s = ""
    for i in 0..<m.stride:
      s.addFloat m.p.data[e + i]
      if i < m.stride - 1:
        s.add "\t"
    result.add s

proc main =
  var a = matrix(@[1.0, 2, 3, 4, 5, 6, 7, 8], 4)
  echo a
  echo a[3, 0]

main()
