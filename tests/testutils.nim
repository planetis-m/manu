import math, morpheus / matrix

# private utility routines

proc check*(x, y: float) =
   # Check magnitude of difference of scalars.
   let eps = pow(2.0, -52.0)
   if x == 0.0 and abs(y) < 10.0 * eps: return
   if y == 0.0 and abs(x) < 10.0 * eps: return
   if abs(x-y) > 10.0 * eps * max(abs(x), abs(y)):
      raise newException(ValueError, "The difference x-y is too large: x = " & $x & "  y = " & $y)

proc check*(x, y: seq[float]) =
   # Check norm of difference of "vectors".
   if x.len == y.len:
      for i in 0 ..< x.len:
         check(x[i], y[i])
   else:
      raise newException(ValueError, "Attempt to compare vectors of different lengths")

proc check*(a, b: Matrix) =
   # Check norm of difference of Matrices.
   let eps = pow(2.0, -52.0)
   let x_norm1 = a.norm1()
   let y_norm1 = b.norm1()
   let xmiy_norm1 = norm1(a - b)
   if x_norm1 == 0.0 and y_norm1 < 10.0 * eps: return
   if y_norm1 == 0.0 and x_norm1 < 10.0 * eps: return
   if xmiy_norm1 > 1000.0 * eps * max(x_norm1, y_norm1):
      raise newException(ValueError, "The norm of (a-b) is too large: " & $xmiy_norm1)

proc check*(x, y: seq[seq[float]]) =
   # Check norm of difference of arrays.
   let a = matrix(x)
   let b = matrix(y)
   check(a, b)

proc try_success*(s, e: string) =
   # Print appropriate messages for successful outcome try
   echo(">    ", s, "success")
   if e != "":
      echo(">      Message: ", e)

proc try_failure*(count: var int, s, e: string) =
   # Print appropriate messages for unsuccessful outcome try
   echo(">    ", s, "*** failure ***\n>      Message: ", e)
   inc(count)

proc try_warning*(count: var int, s, e: string) =
   # Print appropriate messages for unsuccessful outcome try
   echo(">    ", s, "*** warning ***\n>      Message: ", e)
   inc(count)
