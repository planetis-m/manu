type
   ColVector* = distinct Matrix
   RowVector* = distinct Matrix

template colVector*(m: int): Matrix = ColVector(matrix(m, 1))
template rowVector*(n: int): Matrix = RowVector(matrix(1, n))

proc `[]`*(v: ColVector, i: int): float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < v.m)
   m.data[i]

proc `[]`*(v: var ColVector, i: int): var float {.inline.} =
   ## Get a single element.
   checkBounds(i >= 0 and i < v.m)
   m.data[i]

proc `[]=`*(v: var ColVector, i: int, s: float) {.inline.} =
   ## Set a single element.
   checkBounds(i >= 0 and i < v.m)
   m.data[i] = s

proc `[]`*(v: RowVector, j: int): float {.inline.} =
   ## Get a single element.
   checkBounds(j >= 0 and j < v.n)
   m.data[j]

proc `[]`*(v: var RowVector, j: int): var float {.inline.} =
   ## Get a single element.
   checkBounds(j >= 0 and j < v.n)
   m.data[j]

proc `[]=`*(v: var RowVector, j: int, s: float) {.inline.} =
   ## Set a single element.
   checkBounds(j >= 0 and j < v.n)
   m.data[j] = s

proc `+`*(a: sink Matrix, b: ColVector): Matrix =
   ## ``C = A + B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[i]

proc `+`*(a: sink Matrix, b: RowVector): Matrix =
   ## ``C = A + B``
   assert(b.n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] + b[j]

proc `+`*(a: ColVector, b: RowVector): Matrix =
   ## ``C = A + B``
   result = matrix(a.m, b.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] + b[j]

template `+`*(b: RowVector, a: ColVector): Matrix = a + b
template `+`*(b: ColVector, a: Matrix): Matrix = m + b
template `+`*(b: RowVector, a: Matrix): Matrix = m + b

proc `+=`*(a: var Matrix, b: ColVector) =
   ## ``A = A + B``
   assert(b.m == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[i]

proc `+=`*(a: var Matrix, b: RowVector) =
   ## ``A = A + B``
   assert(b.n == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] + b[i]

proc `-`*(a: sink Matrix, b: ColVector): Matrix =
   ## ``C = A - B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[i]

proc `-`*(a: sink Matrix, b: RowVector): Matrix =
   ## ``C = A - B``
   assert(b.n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] - b[j]

proc `-`*(a: ColVector, b: RowVector): Matrix =
   ## ``C = A - B``
   result = matrix(a.m, b.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] - b[j]

template `-`*(b: RowVector, a: ColVector): Matrix = a - b
template `-`*(b: ColVector, a: Matrix): Matrix = m - b
template `-`*(b: RowVector, a: Matrix): Matrix = m - b

proc `-=`*(a: var Matrix, b: ColVector) =
   ## ``A = A - B``
   assert(b.m == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[i]

proc `-=`*(a: var Matrix, b: RowVector) =
   ## ``A = A - B``
   assert(b.n == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] - b[i]

proc `*.`*(a: sink Matrix, b: ColVector): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[i]

proc `*.`*(a: sink Matrix, b: RowVector): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   assert(b.n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] * b[j]

proc `*.`*(a: ColVector, b: RowVector): Matrix =
   ## Element-by-element multiplication, ``C = A.*B``
   result = matrix(a.m, b.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] * b[j]

template `*.`*(b: RowVector, a: ColVector): Matrix = a *. b
template `*.`*(b: ColVector, a: Matrix): Matrix = m *. b
template `*.`*(b: RowVector, a: Matrix): Matrix = m *. b

proc `*.=`*(a: var Matrix, b: ColVector) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.m == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[i]

proc `*.=`*(a: var Matrix, b: RowVector) =
   ## Element-by-element multiplication in place, ``A = A.*B``
   assert(b.n == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] * b[i]

proc `/.`*(a: sink Matrix, b: ColVector): Matrix =
   ## Element-by-element right division, ``C = A./B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / b[i]

proc `/.`*(a: sink Matrix, b: RowVector): Matrix =
   ## Element-by-element right division, ``C = A./B``
   assert(b.n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = result[i, j] / b[j]

proc `/.`*(a: ColVector, b: RowVector): Matrix =
   ## Element-by-element right division, ``C = A./B``
   result = matrix(a.m, b.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = a[i] / b[j]

proc `/.=`*(a: var Matrix, b: ColVector) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.m == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[i]

proc `/.=`*(a: var Matrix, b: RowVector) =
   ## Element-by-element right division in place, ``A = A./B``
   assert(b.n == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = a[i, j] / b[i]

proc `\.`*(a: sink Matrix, b: ColVector): Matrix =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.m == a.m, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[i] / result[i, j]

proc `\.`*(a: sink Matrix, b: RowVector): Matrix =
   ## Element-by-element left division, ``C = A.\B``
   assert(b.n == a.n, "Matrix-vector dimensions must agree.")
   result = a
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[j] / result[i, j]

proc `\.`*(a: ColVector, b: RowVector): Matrix =
   ## Element-by-element left division, ``C = A.\B``
   result = matrix(a.m, b.n)
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = b[j] / a[i]

proc `\.=`*(a: var Matrix, b: ColVector) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.m == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[i] / a[i, j]

proc `\.=`*(a: var Matrix, b: RowVector) =
   ## Element-by-element left division in place, ``A = A.\B``
   assert(b.n == a.m, "Matrix dimensions must agree.")
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = b[i] / a[i, j]

template `/.`*(b: RowVector, a: ColVector): Matrix = a \. b
template `/.`*(b: ColVector, a: Matrix): Matrix = m \. b
template `/.`*(b: RowVector, a: Matrix): Matrix = m \. b

template `\.`*(b: RowVector, a: ColVector): Matrix = a /. b
template `\.`*(b: ColVector, a: Matrix): Matrix = m /. b
template `\.`*(b: RowVector, a: Matrix): Matrix = m /. b

proc sumColumns*(m: Matrix): RowVector =
   ## Column sum.
   ##
   ## ``return``: An 1-by-n matrix with the column sum of each row.
   result = rowVector(m.n)
   for j in 0 ..< m.n:
      var s = 0.0
      for i in 0 ..< m.m:
         s += m[i, j]
      result[j] = s

proc sumRows*(m: Matrix): ColVector =
   ## Row sum.
   ##
   ## ``return``: An m-by-1 matrix with the row sum of each column.
   result = colVector(m.m)
   for i in 0 ..< m.m:
      var s = 0.0
      for j in 0 ..< m.n:
         s += m[i, j]
      result[i] = s

template `$`*(v: ColVector): string = $Matrix(v)
template `$`*(v: RowVector): string = $Matrix(v)
