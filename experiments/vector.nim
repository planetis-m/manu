type
   Vector* = object
      m: int
      data: seq[float] # Array for internal storage of elements.

proc vector*(m: int): Vector {.inline.} =
   ## Construct an m vector of zeros.
   result.m = m
   result.data = newSeq[float](m)

proc vector*(m: int, s: float): Vector =
   ## Construct an m constant vector.
   result.m = m
   result.data = newSeqUninitialized[float](m)
   for i in 0 ..< result.data.len:
      result.data[i] = s

template ones*(m: int): Matrix = vector(m, 1.0)
template zeros*(m: int): Matrix = vector(m)

proc vector*(data: sink seq[float]): Vector =
   ## Construct a vector from a one-dimensional array.
   ##
   ## parameter ``data``: one-dimensional array of float.
   result.m = data.len
   result.data = data

proc randVector*(m: int, max: float): Vector =
   ## Generate vector with random elements.
   ##
   ## ``return``: an m vector with uniformly distributed random elements.
   result.m = m
   result.data = newSeqUninitialized[float](result.m)
   for i in 0 ..< result.m:
      result.data[i] = rand(max)

proc randVector*(m: int, x: Slice[float]): Vector =
   ## Generate vector with random elements.
   ##
   ## ``return``: an m vector with uniformly distributed random elements.
   result.m = m
   result.data = newSeqUninitialized[float](result.m)
   for i in 0 ..< result.m:
      result.data[i] = rand(x)

template randVector*(m, n: int): Vector = randVector(m, 1.0)

proc getArray*(v: Vector): seq[float] {.inline.} =
   ## Copy the internal one-dimensional array.
   v.data

proc m*(v: Vector): int {.inline.} =
   ## Get dimension.
   v.m

proc `[]`*(v: Vector, i: int): float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < v.m)
   v.data[i]

proc `[]`*(v: var Vector, i: int): var float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < v.m)
   v.data[i]

proc `[]=`*(v: var Vector, i: int, s: float) {.inline.} =
   ## Set a single element.
   checkBounds(i >= 0 and i < v.m)
   v.data[i] = s

template `^^`(dim, i: untyped): untyped =
  (when i is BackwardsIndex: dim - int(i) else: int(i))

proc `[]`*[U, V](v: Vector, r: HSlice[U, V]): Vector =
   ## Get a subvector,
   ## ``v[i0 .. i1]``
   let ra = v.m ^^ r.a
   let rb = v.m ^^ r.b
   checkBounds(ra >= 0 and rb < v.m, "Subvector dimensions")
   result = vector(rb - ra + 1)
   for i in 0 ..< result.m:
      result[i] = v[i + ra]

proc `[]`*(v: Vector, r: openarray[int]): Vector =
   ## Get a subvector,
   ## ``v[[0, 2, 3, 4]]``
   checkBounds(r.len <= v.m, "Subvector dimensions")
   result = vector(r.len)
   for i in 0 ..< r.len:
      result[i] = v[r[i]]

proc `[]=`*[U, V](v: var Vector, r: HSlice[U, V], a: Vector) =
   ## Set a subvector,
   ## ``m[i0 .. i1] = a``
   let ra = v.m ^^ r.a
   let rb = v.m ^^ r.b
   checkBounds(rb - ra + 1 == a.m, "Subvector dimensions")
   for i in 0 ..< a.m:
      v[i + ra] = a[i]

proc `[]=`*(v: var Vector, r: openarray[int], a: Vector) =
   ## Set a subvector,
   ## ``m[[0, 2, 3, 4]] = a``
   checkBounds(r.len == a.m, "Subvector dimensions")
   for i in 0 ..< r.len:
      v[r[i]] = a[i]

proc `-`*(v: sink Vector): Vector =
   ## Unary minus
   result = v
   for i in 0 ..< result.m:
      result[i] = -result[i]

proc `+`*(a: sink Vector, b: Vector): Vector =
   ## ``C = A + B``
   assert(b.m == a.m , "Vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      result[i] = result[i] + b[i]

proc `+=`*(a: var Vector, b: Vector) =
   ## ``A = A + B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   for i in 0 ..< a.m:
      a[i] = a[i] + b[i]

proc `-`*(a: sink Vector, b: Vector): Vector =
   ## ``C = A - B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      result[i] = result[i] - b[i]

proc `-=`*(a: var Vector, b: Vector) =
   ## ``A = A - B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   for i in 0 ..< a.m:
      a[i] = a[i] - b[i]

proc `*`*(a: sink Vector, b: Vector): Vector =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      result[i] = result[i] * b[i]

proc `*=`*(a: var Vector, b: Vector) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   for i in 0 ..< a.m:
      a[i] = a[i] * b[i]

proc `/`*(a: sink Vector, b: Vector): Vector =
   ## Element-by-element right division, ``C = A./B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      result[i] = result[i] / b[i]

proc `/=`*(a: var Vector, b: Vector) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   for i in 0 ..< a.m:
      a[i] = a[i] / b[i]

proc `\`*(a: sink Vector, b: Vector): Vector =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      result[i] = b[i] / result[i]

proc `\=`*(a: var Vector, b: Vector) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.m == a.m, "Vector dimensions must agree.")
   for i in 0 ..< a.m:
      a[i] = b[i] / a[i]

proc `+`*(v: sink Vector, s: float): Vector =
   ## Add a vector to a scalar, ``C = s+A``
   result = v
   for i in 0 ..< result.m:
      result[i] = s + result[i]

template `+`*(s: float, v: Vector): Vector = v + s
template `-`*(v: Vector, s: float): Vector = v + (-s)

proc `+=`*(v: var Vector, s: float) =
   ## Add a vector to a scalar in place, ``A = s+A``
   for i in 0 ..< v.m:
      v[i] = s + v[i]

template `-=`*(v: Vector, s: float): Vector = m += (-s)

proc `-`*(s: float, v: sink Vector): Vector =
   ## Subtract a vector from a scalar, ``C = s-A``
   result = v
   for i in 0 ..< result.m:
      result[i] = s - result[i]

proc `*`*(v: sink Vector, s: float): Vector =
   ## Multiply a matrix by a scalar, ``C = s*A``
   result = v
   for i in 0 ..< result.m:
      result[i] = s * result[i]

proc `*=`*(v: var Vector, s: float) =
   ## Multiply a vector by a scalar in place, ``A = s*A``
   for i in 0 ..< v.m:
      v[i] = s * v[i]

template `*`*(s: float, v: Vector): Vector = v * s

proc `/`*(v: sink Vector, s: float): Vector =
   ## Divide a vector by a scalar, ``C = A/s``
   result = v
   for i in 0 ..< result.m:
      result[i] = result[i] / s

template `/`*(s: float, v: Vector): Vector = v * (1 / s)

proc `/=`*(v: var Vector, s: float) =
   ## Divide a matrix by a scalar in place, ``A = A/s``
   for i in 0 ..< v.m:
      v[i] = v[i] / s

proc sum*(v: Vector): float =
   ## Sum.
   for i in 0 ..< v.m:
      result += v[i]

proc `+-`*(a: sink Matrix, b: Vector): Matrix =
   ## ``C = A + B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[i]

proc `+|`*(a: sink Matrix, b: Vector): Matrix =
   ## ``C = A + B``
   assert(b.m == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[j]

proc `--`*(a: sink Matrix, b: Vector): Matrix = # this is ridiculous
   ## ``C = A + B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[i]

proc `-|`*(a: sink Matrix, b: Vector): Matrix =
   ## ``C = A + B``
   assert(b.m == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[j]

proc sumColumns*(m: Matrix): Vector =
   ## Column sum.
   ##
   ## ``return``: An n vector with the column sum of each row.
   result = vector(m.n)
   for j in 0 ..< m.n:
      var s = 0.0
      for i in 0 ..< m.m:
         s += m[i, j]
      result[j] = s

proc sumRows*(m: Matrix): Vector =
   ## Row sum.
   ##
   ## ``return``: An m vector with the row sum of each column.
   result = vector(m.m)
   for i in 0 ..< m.m:
      var s = 0.0
      for j in 0 ..< m.n:
         s += m[i, j]
      result[i] = s

proc `$`*(v: Vector): string =
   var formatted = newSeqOfCap[string](v.m)
   var mColj = newSeq[float](v.m)
   for i in 0 ..< v.m:
      mColj[i] = v[i]
   formatted.add columnFormat(mColj)
   result = ""
   for i in 0 ..< v.m:
      if i == 0:
         result.add "⎡"
      elif i == v.m - 1:
         result.add "⎣"
      else:
         result.add "⎢"
      result.add formatted[i]
      if i == 0:
         result.add "⎤\n"
      elif i == v.m - 1:
         result.add "⎦"
      else:
         result.add "⎥\n"

template makeUniversal*(fname: untyped) =
   proc fname*(v: sink Vector): Vector =
      result = v
      for i in 0 ..< v.m:
         result.data[i] = fname(result.data[i])
   proc fname*(m: sink Matrix): Matrix =
      let len = m.m * m.n
      result = m
      for i in 0 ..< len:
         result.data[i] = fname(result.data[i])

template makeUniversalBinary*(fname: untyped) =
   proc fname*(a: sink Vector, b: Vector): Vector =
      assert(b.m == a.m, "Vector dimensions must agree.")
      result = a
      for i in 0 ..< result.m:
         result.data[i] = fname(result.data[i], b.data[i])
   proc fname*(a: sink Matrix, b: Matrix): Matrix =
      assert(b.m == a.m and b.n == a.n, "Matrix dimensions must agree.")
      let len = a.m * a.n
      result = a
      for i in 0 ..< len:
         result.data[i] = fname(result.data[i], b.data[i])
