import manu, utils

# main
proc main() =
  var A, B, C, Z, O, R, S: Matrix[float32]
  var errorCount = 0
  var warningCount = 0
  let columnwise = @[1'f32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
  let nonconformld = 4 # leading dimension which is valid, but nonconforming
  let validld = 3 # leading dimension of intended test Matrices

  A = matrix(columnwise, validld)
  # Array-like methods:
  echo("\nTesting array-like methods...")
  S = matrix(columnwise, nonconformld)
  R = randMatrix[float32](A.rowDimension, A.columnDimension)
  A = R
  try:
    S = A - S
    tryFailure(errorCount, "minus conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("minus conformance check... ", "")
  if (A - R).norm1() != 0:
    tryFailure(errorCount, "minus... ",
                "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)")
  else:
    trySuccess("minus... ", "")
  A = R #.copy()
  A -= R
  Z = matrix[float32](A.rowDimension, A.columnDimension)
  try:
    A -= S
    tryFailure(errorCount, "minusEquals conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("minusEquals conformance check... ", "")
  if (A - Z).norm1() != 0:
    tryFailure(errorCount, "minusEquals... ",
                "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)")
  else:
    trySuccess("minusEquals... ", "")
  A = R #.copy()
  B = randMatrix[float32](A.rowDimension, A.columnDimension)
  C = A - B
  try:
    S = A + S
    tryFailure(errorCount, "plus conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("plus conformance check... ", "")
  try:
    check(C + B, A)
    trySuccess("plus... ", "")
  except ValueError:
    tryFailure(errorCount, "plus... ", "(C = A - B, but C + B != A)")
  C = A - B
  C += B
  try:
    A += S
    tryFailure(errorCount, "plusEquals conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("plusEquals conformance check... ", "")
  try:
    check(C, A)
    trySuccess("plusEquals... ", "")
  except ValueError:
    tryFailure(errorCount, "plusEquals... ", "(C = A - B, but C = C + B != A)")
  A = -R
  try:
    check(A + R, Z)
    trySuccess("uminus... ", "")
  except ValueError:
    tryFailure(errorCount, "uminus... ", "(-A + A != zeros)")
  A = R #.copy()
  O = matrix(A.rowDimension, A.columnDimension, 1'f32)
  #C = A \. R
  #try:
    #S = A \. S
    #tryFailure(errorCount, "arrayLeftDivide conformance check... ", "nonconformance not raised")
  #except AssertionDefect:
    #trySuccess("arrayLeftDivide conformance check... ", "")
  #try:
    #check(C, O)
    #trySuccess("arrayLeftDivide... ", "")
  #except ValueError:
    #tryFailure(errorCount, "arrayLeftDivide... ", "(M.\\M != ones)")
  #try:
    #A \.= S
    #tryFailure(errorCount, "arrayLeftDivideEquals conformance check... ", "nonconformance not raised")
  #except AssertionDefect:
    #trySuccess("arrayLeftDivideEquals conformance check... ", "")
  #A \.= R
  #try:
    #check(A,O)
    #trySuccess("arrayLeftDivideEquals... ", "")
  #except ValueError:
    #tryFailure(errorCount, "arrayLeftDivideEquals... ", "(M.\\M != ones)")
  #A = R #.copy()
  try:
    C = A /. S
    tryFailure(errorCount, "arrayRightDivide conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("arrayRightDivide conformance check... ", "")
  C = A /. R
  try:
    check(C, O)
    trySuccess("arrayRightDivide... ", "")
  except ValueError:
    tryFailure(errorCount, "arrayRightDivide... ", "(M./M != ones)")
  try:
    A /.= S
    tryFailure(errorCount, "arrayRightDivideEquals conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("arrayRightDivideEquals conformance check... ", "")
  A /.= R
  try:
    check(A, O)
    trySuccess("arrayRightDivideEquals... ", "")
  except ValueError:
    tryFailure(errorCount, "arrayRightDivideEquals... ", "(M./M != ones)")
  A = R #.copy()
  B = randMatrix[float32](A.rowDimension, A.columnDimension)
  try:
    S = A *. S
    tryFailure(errorCount, "arrayTimes conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("arrayTimes conformance check... ", "")
  C = A *. B
  try:
    check(C /. B, A)
    trySuccess("arrayTimes... ", "")
  except ValueError:
    tryFailure(errorCount, "arrayTimes... ", "(A = R, C = A.*B, but C./B != A)")
  try:
    A *.= S
    tryFailure(errorCount, "arrayTimesEquals conformance check... ", "nonconformance not raised")
  except AssertionDefect:
    trySuccess("arrayTimesEquals conformance check... ", "")
  A *.= B
  try:
    check(A /. B, R)
    trySuccess("arrayTimesEquals... ", "")
  except ValueError:
    tryFailure(errorCount, "arrayTimesEquals... ", "(A = R, A = A.*B, but A./B != R)")

  echo("\nTestMatrixOps completed.")
  echo("Total errors reported: ", errorCount)
  #echo("Total warnings reported: ", warningCount)
  if errorCount > 0: raise newException(ValueError, "")

try:
  main()
except:
  quit(QuitFailure)
