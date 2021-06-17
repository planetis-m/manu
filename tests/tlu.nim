import utils, manu / [matrix, lu]

proc main =
  var errorCount = 0
  var warningCount = 0

  echo("Testing LU Decomposition")
  block:
    let
      columnwise = @[14.0, 7, 6, 18]
      mat = matrix(columnwise, 2)
      lu = mat.lu()
    try:
      check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
      trySuccess("LUDecomposition 2x2...", "")
    except ValueError:
      tryFailure(errorCount, "LUDecomposition 2x2...",
          "incorrect LU decomposition calculation")

  block:
    let
      columnwise = @[1.0, 0, 2, 0, 10, 0, 2, 0, 9]
      mat = matrix(columnwise, 3)
      lu = mat.lu()
    try:
      check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
      trySuccess("LUDecomposition 3x3...", "")
    except ValueError:
      tryFailure(errorCount, "LUDecomposition 3x3...",
          "incorrect LU decomposition calculation")

  block:
    let
      columnwise = @[77.0, 19, 1, -9, 17, 34, 11, 100, -2]
      mat = matrix(columnwise, 3)
      lu = mat.lu()
    try:
      check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
      trySuccess("LUDecomposition 3x3_2...", "")
    except ValueError:
      tryFailure(errorCount, "LUDecomposition 3x3_2...",
          "incorrect LU decomposition calculation")

  block:
    let
      columnwise = @[99.0, 14, 39, 11, 1, 65, 40, 5, -10, 7, -2, 43, 6, 48, 9, 99]
      mat = matrix(columnwise, 4)
      lu = mat.lu()
    try:
      check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
      trySuccess("LUDecomposition 4x4...", "")
    except ValueError:
      tryFailure(errorCount, "LUDecomposition 4x4...",
          "incorrect LU decomposition calculation")

  block:
    let
      columnwise = @[6.0, 19, 81, 10, 65, 100, 1, -10, 16, 71, 58, -17, 88, 19,
                      29, -44, 4, 16, -100, 1, 50, 76, 93, 35, -24]
      mat = matrix(columnwise, 5)
      lu = mat.lu()
    try:
      check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
      trySuccess("LUDecomposition 5x5...", "")
    except ValueError:
      tryFailure(errorCount, "LUDecomposition 5x5...",
          "incorrect LU decomposition calculation")

  block:
    let
      columnwise = @[88.0, 17, 6, 14, -1, 5, -5, 41, 16, -29, 7, -53, 19, 22,
                    -99, 3, 101, -91, 8, 26, 71, -66, 46, 18, 23]
      mat = matrix(columnwise, 5)
      lu = mat.lu()
    try:
      check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
      trySuccess("LUDecomposition 5x5_2...", "")
    except ValueError:
      tryFailure(errorCount, "LUDecomposition 5x5_2...",
          "incorrect LU decomposition calculation")

  echo("\nTestLU completed.")
  echo("Total errors reported: ", errorCount)
  echo("Total warnings reported: ", warningCount)
  if errorCount > 0: raise newException(ValueError, "")

try:
  main()
except:
  quit(QuitFailure)
