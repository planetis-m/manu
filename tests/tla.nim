## TestMatrix tests the functionality of the Jama Matrix class and associated decompositions.
##
## Detailed output is provided indicating the functionality being tested
## and whether the functionality is correctly implemented.   Exception handling
## is also tested.
##
## The test is designed to run to completion and give a summary of any implementation errors
## encountered. The final output should be:
## .. code-block::
##   TestMatrix completed.
##   Total errors reported: n1
##   Total warning reported: n2
##
## If the test does not run to completion, this indicates that there is a
## substantial problem within the implementation that was not anticipated in the test design.
## The stopping point should give an indication of where the problem exists.
import std/math, manu, utils

# main
proc main() =
  var A, B, Z, O, R, X, SUB, T, SQ, DEF, SOL: Matrix[float32]
  var errorCount = 0
  var warningCount = 0
  let columnwise = @[1'f32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
  let avals = @[@[1'f32, 4, 7, 10], @[2'f32, 5, 8, 11], @[3'f32, 6, 9, 12]]
  let subavals = @[@[5'f32, 8, 11], @[6'f32, 9, 12]]
  let rankdef = avals
  let tvals = @[@[1'f32, 2, 3], @[4'f32, 5, 6], @[7'f32, 8, 9], @[10'f32, 11, 12]]
  let pvals = @[@[4'f32, 1, 1], @[1'f32, 2, 3], @[1'f32, 3, 6]]
  let evals =
    @[@[0'f32, 1, 0, 0], @[1'f32, 0, 2e-7, 0], @[0'f32, -2e-7, 0, 1], @[0'f32, 0, 1, 0]]
  let square = @[@[166'f32, 188, 210], @[188'f32, 214, 240], @[210'f32, 240, 270]]
  let sqSolution = @[@[13'f32], @[15'f32]]
  let condmat = @[@[1'f32, 3], @[7'f32, 9]]
  let badeigs = @[@[0'f32, 0, 0, 0, 0], @[0'f32, 0, 0, 0, 1], @[0'f32, 0, 0, 1, 0],
     @[1'f32, 1, 0, 0, 1], @[1'f32, 0, 1, 0, 1]]
  let columnsummax = 33'f32
  let rowsummax = 30'f32
  let sumofdiagonals = 15'f32
  let sumofsquares = 650'f32

  # LA methods:
  #   transpose
  #   times
  #   cond
  #   rank
  #   det
  #   trace
  #   norm1
  #   norm2
  #   normF
  #   normInf
  #   solve
  #   solveTranspose
  #   inverse
  #   chol
  #   eig
  #   lu
  #   qr
  #   svd
  echo("\nTesting linear algebra methods...")
  SUB = matrix(subavals)
  A = matrix(columnwise, 3)
  Z = matrix[float32](A.rowDimension, A.columnDimension)
  T = matrix(tvals)
  T = A.transpose()
  try:
    check(A.transpose(), T)
    trySuccess("transpose...", "")
  except ValueError:
    tryFailure(errorCount, "transpose()...", "transpose unsuccessful")
  try:
    check(A.norm1(), columnsummax)
    trySuccess("norm1...", "")
  except ValueError:
    tryFailure(errorCount, "norm1()...", "incorrect norm calculation")
  try:
    check(A.normInf(), rowsummax)
    trySuccess("normInf()...", "")
  except ValueError:
    tryFailure(errorCount, "normInf()...", "incorrect norm calculation")
  try:
    check(A.normF(), sqrt(sumofsquares))
    trySuccess("normF...", "")
  except ValueError:
    tryFailure(errorCount, "normF()...", "incorrect norm calculation")
  try:
    check(A.trace(), sumofdiagonals)
    trySuccess("trace()...", "")
  except ValueError:
    tryFailure(errorCount, "trace()...", "incorrect trace calculation")
  try:
    check(A[0 .. A.rowDimension-1, 0 .. A.rowDimension-1].det(), 0)
    trySuccess("det()...", "")
  except ValueError:
    tryFailure(errorCount, "det()...",
                "incorrect determinant calculation")
  SQ = matrix(square)
  try:
    check(A * A.transpose(), SQ)
    trySuccess("times(Matrix)...", "")
  except ValueError:
    tryFailure(errorCount, "times(Matrix)...",
                "incorrect Matrix-Matrix product calculation")
  try:
    check(A * 0, Z)
    trySuccess("times(double)...", "")
  except ValueError:
    tryFailure(errorCount, "times(double)...",
                "incorrect Matrix-scalar product calculation")
  A = matrix(columnwise, 4)
  let QR = A.qr()
  R = QR.getR()
  try:
    check(A, QR.getQ() * R)
    trySuccess("QRDecomposition...", "")
  except ValueError:
    tryFailure(errorCount, "QRDecomposition...", "incorrect QR decomposition calculation")
  var SVD = A.svd()
  try:
    check(A, SVD.getU() * SVD.getS() * SVD.getV().transpose())
    trySuccess("SingularValueDecomposition...", "")
  except ValueError:
    tryFailure(errorCount, "SingularValueDecomposition...",
                "incorrect singular value decomposition calculation")
  DEF = matrix(rankdef)
  try:
    assert(DEF.rank() == min(DEF.rowDimension, DEF.columnDimension)-1)
    trySuccess("rank()...", "")
  except AssertionDefect:
    tryFailure(errorCount, "rank()...", "incorrect rank calculation")
  B = matrix(condmat)
  SVD = B.svd()
  let singularvalues = SVD.getSingularValues()
  try:
    check(B.cond(), singularvalues[0] / singularvalues[min(B.rowDimension,
          B.columnDimension)-1])
    trySuccess("cond()...", "")
  except ValueError:
    tryFailure(errorCount, "cond()...", "incorrect condition number calculation")
  let n = A.columnDimension
  A = A[0 .. n-1, 0 .. n-1]
  A[0, 0] = 0
  let LU = A.lu()
  try:
    check(A[LU.getPivot(), 0 .. n-1], LU.getL() * LU.getU())
    trySuccess("LUDecomposition...", "")
  except ValueError:
    tryFailure(errorCount, "LUDecomposition...", "incorrect LU decomposition calculation")
  X = A.inverse()
  try:
    check(A * X, identity[float32](3))
    trySuccess("inverse()...", "")
  except ValueError:
    tryFailure(errorCount, "inverse()...", "incorrect inverse calculation")
  O = matrix(SUB.rowDimension, 1, 1'f32)
  SOL = matrix(sqSolution)
  SQ = SUB[0 .. SUB.rowDimension-1, 0 .. SUB.rowDimension-1]
  try:
    check(SQ.solve(SOL), O)
    trySuccess("solve()...", "")
  except AssertionDefect as e1:
    tryFailure(errorCount, "solve()...", e1.msg)
  except ValueError as e:
    tryFailure(errorCount, "solve()...", e.msg)
  A = matrix(pvals)
  let Chol = A.chol()
  let L = Chol.getL()
  try:
    check(A, L * L.transpose())
    trySuccess("CholeskyDecomposition...", "")
  except ValueError:
    tryFailure(errorCount, "CholeskyDecomposition...",
                "incorrect Cholesky decomposition calculation")
  X = Chol.solve(identity[float32](3))
  try:
    check(A * X, identity[float32](3))
    trySuccess("CholeskyDecomposition solve()...", "")
  except ValueError:
    tryFailure(errorCount, "CholeskyDecomposition solve()...",
                "incorrect Choleskydecomposition solve calculation")
  var Eig = A.eig()
  var D = Eig.getD()
  var V = Eig.getV()
  try:
    check(A * V, V * D)
    trySuccess("EigenvalueDecomposition (symmetric)...", "")
  except ValueError:
    tryFailure(errorCount, "EigenvalueDecomposition (symmetric)...",
                "incorrect symmetric Eigenvalue decomposition calculation")
  A = matrix(evals)
  Eig = A.eig()
  D = Eig.getD()
  V = Eig.getV()
  try:
    check(A * V, V * D)
    trySuccess("EigenvalueDecomposition (nonsymmetric)...", "")
  except ValueError:
    tryFailure(errorCount, "EigenvalueDecomposition (nonsymmetric)...",
                "incorrect nonsymmetric Eigenvalue decomposition calculation")
  try:
    echo("\nTesting Eigenvalue; If this hangs, we've failed")
    let bA = matrix(badeigs)
    let bEig = bA.eig()
    trySuccess("EigenvalueDecomposition (hang)...", "")
  except ValueError:
    tryFailure(errorCount, "EigenvalueDecomposition (hang)...",
                "incorrect termination")

  echo("\nTestLA completed.")
  echo("Total errors reported: ", errorCount)
  #echo("Total warnings reported: ", warningCount)
  if errorCount > 0: raise newException(ValueError, "")

try:
  main()
except:
  quit(QuitFailure)
