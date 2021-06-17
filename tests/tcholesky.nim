import utils, manu / [matrix, cholesky]

proc main =
  var errorCount = 0
  var warningCount = 0

  echo("Testing Cholesky Decomposition")
  block:
    let
      columnwise = @[1.0, -2, -2, 5]
      mat = matrix(columnwise, 2)
      l = chol(mat).getL
    try:
      check(l * l.transpose, mat)
      trySuccess("CholeskyDecomposition 2x2...", "")
    except ValueError:
      tryFailure(errorCount, "CholeskyDecomposition 2x2...",
          "incorrect Cholesky decomposition calculation")

  block:
    let
      columnwise = @[25.0, 15, -5, 15, 18, 0, -5, 0, 11]
      mat = matrix(columnwise, 3)
      l = chol(mat).getL
    try:
      check(l * l.transpose, mat)
      trySuccess("CholeskyDecomposition 3x3...", "")
    except ValueError:
      tryFailure(errorCount, "CholeskyDecomposition 3x3...",
          "incorrect Cholesky decomposition calculation")

  block:
    let
      columnwise = @[18.0, 22, 54, 42, 22, 70, 86, 62, 54, 86, 174, 134, 42, 62,
          134, 106]
      mat = matrix(columnwise, 4)
      l = chol(mat).getL
    try:
      check(l * l.transpose, mat)
      trySuccess("CholeskyDecomposition 4x4...", "")
    except ValueError:
      tryFailure(errorCount, "CholeskyDecomposition 4x4...",
          "incorrect Cholesky decomposition calculation")

  block:
    let
      columnwise = @[6.0, 15, 55, 72, 101, 15, 55, 225, 229, 256, 55, 225, 979,
                    1024, 1200, 72, 229, 1024, 2048, 2057, 101, 256, 1200, 2057, 6000]
      mat = matrix(columnwise, 5)
      l = chol(mat).getL
    try:
      check(l * l.transpose, mat)
      trySuccess("CholeskyDecomposition 5x5...", "")
    except ValueError:
      tryFailure(errorCount, "CholeskyDecomposition 5x5...",
          "incorrect Cholesky decomposition calculation")

  echo("\nTestCholesky completed.")
  echo("Total errors reported: ", errorCount)
  echo("Total warnings reported: ", warningCount)
  if errorCount > 0: raise newException(ValueError, "")

try:
  main()
except:
  quit(QuitFailure)
