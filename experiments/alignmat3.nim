import memallocs

const
  MemAlign {.define: "manu.memalign".} = 32

template checkBounds(cond: untyped, msg = "") =
  when compileOption("boundChecks"):
    {.line.}:
      if not cond:
        raise newException(IndexDefect, msg)

type
  Matrix*[T: SomeFloat] = object
    m, n: int                # Row and column dimensions.
    p {.noalias.}: ptr UncheckedArray[T] # Array for internal storage of elements.

template allocs(p, s) =
  let p = cast[ptr UncheckedArray[T]](alignedAlloc(s * sizeof(T), MemAlign))

template allocs0(p, s) =
  let p = cast[ptr UncheckedArray[T]](alignedAlloc0(s * sizeof(T), MemAlign))

proc `=destroy`*[T](a: var Matrix[T]) =
  if a.p != nil:
    alignedDealloc(a.p, MemAlign)

template dups(a, b) =
  a.m = b.m
  a.n = b.n
  if b.p != nil:
    let len = b.m * b.n
    allocs(p, len)
    a.p = p
    copyMem(a.p, b.p, len * sizeof(T))

proc `=copy`*[T](a: var Matrix[T]; b: Matrix[T]) =
  if a.p != b.p:
    `=destroy`(a)
    wasMoved(a)
    dups(a, b)

when defined(nimHasDup):
  proc `=dup`*[T](b: Matrix[T]): Matrix[T] =
    dups(result, b)

proc matrix*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
  ## Construct an m-by-n matrix of zeros.
  let len = m * n
  allocs0(p, len)
  result = Matrix[T](m: m, n: n, p: p)

proc `[]`*[T](m: Matrix[T], i, j: int): T {.inline.} =
  ## Get a single element.
  checkBounds(i >= 0 and i < m.m)
  checkBounds(j >= 0 and j < m.n)
  result = m.p[i * m.n + j]

proc `[]`*[T](m: var Matrix[T], i, j: int): var T {.inline.} =
  ## Get a single element.
  checkBounds(i >= 0 and i < m.m)
  checkBounds(j >= 0 and j < m.n)
  m.p[i * m.n + j]

proc `[]=`*[T](m: var Matrix[T], i, j: int, s: T) {.inline.} =
  ## Set a single element.
  checkBounds(i >= 0 and i < m.m)
  checkBounds(j >= 0 and j < m.n)
  m.p[i * m.n + j] = s

proc `*`*[T](a: sink Matrix[T], b: Matrix[T]): Matrix[T] =
  assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
  result = a
  for i in 0 ..< result.m:
    for j in 0 ..< result.n:
      result[i, j] = result[i, j] * b[i, j]

proc main =
  let a, b = matrix[float32](96, 96)
  let c = a * b
  echo c.m

main()
