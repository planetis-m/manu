import manu, utils

proc main() =
  var A, T, Z: Matrix[float32]
  var errorCount = 0
  var warningCount = 0
  let columnwise = @[1'f32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
  A = matrix(columnwise, 3)
  try:
    check(A * 0, Z)
    trySuccess("times(double)...", "")
  except ValueError:
    tryFailure(errorCount, "times(double)...",
                "incorrect Matrix-scalar product calculation")

main()
