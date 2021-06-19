import std/math, manu, utils

proc main() =
  var A, T, Z: Matrix[float32]
  var errorCount = 0
  var warningCount = 0
  let tvals = @[@[1'f32, 2, 3], @[4'f32, 5, 6], @[7'f32, 8, 9], @[10'f32, 11, 12]]
  T = matrix(tvals)
  T = A.transpose()
  try:
    check(A * 0, Z)
    trySuccess("times(double)...", "")
  except ValueError:
    tryFailure(errorCount, "times(double)...",
                "incorrect Matrix-scalar product calculation")

main()
