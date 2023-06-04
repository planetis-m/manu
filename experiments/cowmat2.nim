type
  MatPayloadBase = object
    counter, stride: int

  MatPayload = object
    counter, stride: int
    data: UncheckedArray[float]

  Matrix* = object
    m, n, offset: int
    p: ptr MatPayload

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
  a.offset = b.offset
  a.m = b.m
  a.n = b.n

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
  if y.p == nil:
    result = Matrix(m: 0, n: 0, offset: 0, p: nil)
  else:
    let len = y.m * y.n
    allocs(contentSize(len))
    p.counter = 0
    p.stride = y.n
    # copyMem(addr p.data[0], addr y.p.data[y.offset], len * sizeof(float))
    var rx = 0
    for bx in countup(0, y.p.stride-1):
      for ax in countup(0, y.m*y.p.stride-1, y.p.stride):
        p.data[rx] = y.p.data[y.offset + ax + bx]
        inc rx
    result = Matrix(m: y.m, n: y.n, offset: 0, p: p)

proc matrix*(m, n: int): Matrix {.inline.} =
  ## Construct an m-by-n matrix of zeros.
  let len = m * n
  allocs0(contentSize(len))
  p.stride = n
  result = Matrix(m: m, n: n, offset: 0, p: p)

proc matrix*(data: seq[float], m: int): Matrix =
  ## Construct a matrix from a one-dimensional packed array.
  ##
  ## parameter ``data``: one-dimensional array of float, packed by columns (ala Fortran).
  ## Array length must be a multiple of ``m``.
  let n = if m != 0: data.len div m else: 0
  let len = m * n
  assert(len == data.len, "Array length must be a multiple of m.")
  allocs(contentSize(len))
  p.stride = n
  var rx = 0
  for bx in countup(0, m-1):
    for ax in countup(0, high(data), m):
      p.data[rx] = data[ax + bx]
      inc rx
  result = Matrix(m: m, n: n, offset: 0, p: p)

proc m*(m: Matrix): int {.inline.} =
  ## Get row dimension.
  m.m

proc n*(m: Matrix): int {.inline.} =
  ## Get column dimension.
  m.n

proc `[]`*(m: Matrix, i, j: int): float {.inline.} =
  ## Get a single element.
  checkBounds(i >= 0 and i < m.m)
  checkBounds(j >= 0 and j < m.n)
  result = m.p.data[m.offset + i * m.p.stride + j]

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
  if m.p != nil:
    inc m.p.counter
  result = Matrix(m: rb - ra + 1, n: cb - ca + 1, p: m.p, offset: m.offset + ra * m.n + ca)

proc main =
  let a = matrix(@[1.0, 2, 3, 4, 5, 6, 7, 8], 4)
  let b = a[0..1, 0..1]
  let c = a[1..3, 0..1]
  var d = c[0..2, 1..1]
  d = deepCopy(d)
  var e = c[0..1, 0..0]
  e = deepCopy(e)
  echo (e[0, 0], e[1, 0])
  echo (a[1, 1], b[1, 1], c[0, 1], d[0, 0])

main()
