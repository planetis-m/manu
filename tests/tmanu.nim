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
  var A, B, C, Z, O, I, R, S, X, SUB, M, T, SQ, DEF, SOL: Matrix[float32]
  var errorCount = 0
  var warningCount = 0
  var tmp, s: float32
  let columnwise = @[1'f32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
  let rowwise = @[1'f32, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12]
  let avals = @[@[1'f32, 4, 7, 10], @[2'f32, 5, 8, 11], @[3'f32, 6, 9, 12]]
  let rankdef = avals
  let tvals = @[@[1'f32, 2, 3], @[4'f32, 5, 6], @[7'f32, 8, 9], @[10'f32, 11, 12]]
  let subavals = @[@[5'f32, 8, 11], @[6'f32, 9, 12]]
  var rvals = @[@[2'f32, 5, 8, 11], @[3'f32, 6, 9, 12]]
  rvals.add @[1'f32, 4, 7]
  let pvals = @[@[4'f32, 1, 1], @[1'f32, 2, 3], @[1'f32, 3, 6]]
  let ivals = @[@[1'f32, 0, 0], @[0'f32, 1, 0], @[0'f32, 0, 1]]
  let evals =
    @[@[0'f32, 1, 0, 0], @[1'f32, 0, 2e-7, 0], @[0'f32, -2e-7, 0, 1], @[0'f32,
        0, 1, 0]]
  let square = @[@[166'f32, 188, 210], @[188'f32, 214, 240], @[210'f32, 240, 270]]
  let sqSolution = @[@[13'f32], @[15'f32]]
  let condmat = @[@[1'f32, 3], @[7'f32, 9]]
  let badeigs = @[@[0'f32, 0, 0, 0, 0], @[0'f32, 0, 0, 0, 1], @[0'f32, 0, 0, 1, 0],
     @[1'f32, 1, 0, 0, 1], @[1'f32, 0, 1, 0, 1]]
  let
    rows = 3
    cols = 4
  let invalidld = 5 # should trigger bad shape for construction with val
  let
    raggedr = 0 # (raggedr,raggedc) should be out of bounds in ragged array
    raggedc = 4
  let validld = 3 # leading dimension of intended test Matrices
  let nonconformld = 4 # leading dimension which is valid, but nonconforming
  let
    ib = 1
    ie = 2
    jb = 1
    je = 3 # index ranges for sub Matrix
  let rowindexset = [1, 2]
  let badrowindexset = [1, 3]
  let columnindexset = [1, 2, 3]
  let badcolumnindexset = [1, 2, 4]
  let columnsummax = 33'f32
  let rowsummax = 30'f32
  let sumofdiagonals = 15'f32
  let sumofsquares = 650'f32

  # Constructors
  echo("Testing constructors and constructor-like methods...")
  try:
    # check that exception is thrown in packed constructor with invalid length
    A = matrix(columnwise, invalidld)
    tryFailure(errorCount, "Catch invalid length in packed constructor... ",
                "exception not thrown for invalid input")
  except AssertionDefect as e:
    trySuccess("Catch invalid length in packed constructor... ", e.msg)
  try:
    # check that exception is thrown in default constructor if input array is 'ragged'
    A = matrix(rvals)
    tmp = A[raggedr, raggedc]
  except AssertionDefect as e:
    trySuccess("Catch ragged input to default constructor... ", e.msg)
  except IndexDefect:
    tryFailure(errorCount, "Catch ragged input to constructor... ",
                "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later")
#    try:
#       # check that exception is thrown in constructWithCopy if input array is 'ragged'
#       A = constructWithCopy(rvals)
#       tmp = A[raggedr, raggedc]
#    except AssertionDefect as e:
#       trySuccess("Catch ragged input to constructWithCopy... ", e.msg)
#    except IndexDefect:
#       tryFailure(errorCount, "Catch ragged input to constructWithCopy... ",
#                   "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later")

  A = matrix(columnwise, validld)
  B = matrix(avals)
#    tmp = B[0, 0]
#    avals[0][0] = 0
  C = B - A
#    avals[0][0] = tmp
#    B = constructWithCopy(avals)
#    tmp = B[0, 0]
#    avals[0][0] = 0
#    if tmp - B[0,0] != 0:
#       # check that constructWithCopy behaves properly
#       tryFailure(errorCount, "constructWithCopy... ",
#                   "copy not effected... data visible outside")
#    else:
#       trySuccess("constructWithCopy... ", "")
#
#    avals[0][0] = columnwise[0]
  I = matrix(ivals)
  try:
    check(I, identity[float32](3))
    trySuccess("identity... ", "")
  except ValueError:
    tryFailure(errorCount, "identity... ",
                "identity Matrix not successfully created")

  # Access Methods:
  echo("\nTesting access methods...")

  # Various get methods:
  B = matrix(avals)
  if B.rowDimension != rows:
    tryFailure(errorCount, "rowDimension... ", "")
  else:
    trySuccess("getRowDimension... ", "")

  if B.columnDimension != cols:
    tryFailure(errorCount, "columnDimension... ", "")
  else:
    trySuccess("ColumnDimension... ", "")
  B = matrix(avals)
  var barray = B.getArray()
  if barray != avals:
    tryFailure(errorCount, "getArray... ", "")
  else:
    trySuccess("getArray... ", "")
#    barray = B.getArrayCopy()
#    if barray == avals:
#       tryFailure(errorCount, "getArrayCopy... ", "data not (deep) copied")
#    try:
#       check(barray, avals)
#       trySuccess("getArrayCopy... ", "")
#    except ValueError:
#       tryFailure(errorCount, "getArrayCopy... ", "data not successfully (deep) copied")

  var bpacked = B.getColumnPacked()
  try:
    check(bpacked, columnwise)
    trySuccess("getColumnPacked... ", "")
  except ValueError:
    tryFailure(errorCount, "getColumnPacked... ",
                "data not successfully (deep) copied by columns")
  bpacked = B.getRowPacked()
  try:
    check(bpacked, rowwise)
    trySuccess("getRowPacked... ", "")
  except ValueError:
    tryFailure(errorCount, "getRowPacked... ",
                "data not successfully (deep) copied by rows")
  try:
    tmp = B[B.rowDimension, B.columnDimension-1]
    tryFailure(errorCount, "get(int,int)... ",
                "OutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      tmp = B[B.rowDimension-1, B.columnDimension]
      tryFailure(errorCount, "get(int,int)... ",
                  "OutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("get(int,int)... OutofBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "get(int,int)... ",
                "OutOfBoundsException expected but not thrown")
  try:
    if B[B.rowDimension-1, B.columnDimension-1] !=
          avals[B.rowDimension-1][B.columnDimension-1]:
      tryFailure(errorCount, "get(int,int)... ",
                  "Matrix entry (i,j) not successfully retreived")
    else:
      trySuccess("get(int,int)... ", "")
  except IndexDefect:
    tryFailure(errorCount, "get(int,int)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  SUB = matrix(subavals)
  try:
    M = B[ib .. ie+B.rowDimension+1, jb .. je]
    tryFailure(errorCount, "getMatrix(int,int,int,int)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      M = B[ib .. ie, jb .. je+B.columnDimension+1]
      tryFailure(errorCount, "getMatrix(int,int,int,int)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("getMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "getMatrix(int,int,int,int)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    M = B[ib .. ie, jb .. je]
    try:
      check(SUB, M)
      trySuccess("getMatrix(int,int,int,int)... ", "")
    except ValueError:
      tryFailure(errorCount, "getMatrix(int,int,int,int)... ",
                  "submatrix not successfully retreived")
  except IndexDefect:
    tryFailure(errorCount, "getMatrix(int,int,int,int)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  try:
    M = B[ib .. ie, badcolumnindexset]
    tryFailure(errorCount, "getMatrix(int,int,int[])... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      M = B[ib .. ie+B.rowDimension+1, columnindexset]
      tryFailure(errorCount, "getMatrix(int,int,int[])... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("getMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "getMatrix(int,int,int[])... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    M = B[ib .. ie, columnindexset]
    try:
      check(SUB, M)
      trySuccess("getMatrix(int,int,int[])... ", "")
    except ValueError:
      tryFailure(errorCount, "getMatrix(int,int,int[])... ",
                  "submatrix not successfully retreived")
  except IndexDefect:
    tryFailure(errorCount, "getMatrix(int,int,int[])... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  try:
    M = B[badrowindexset, jb .. je]
    tryFailure(errorCount, "getMatrix(int[],int,int)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      M = B[rowindexset, jb .. je+B.columnDimension+1]
      tryFailure(errorCount, "getMatrix(int[],int,int)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("getMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "getMatrix(int[],int,int)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    M = B[rowindexset, jb .. je]
    try:
      check(SUB, M)
      trySuccess("getMatrix(int[],int,int)... ", "")
    except ValueError:
      tryFailure(errorCount, "getMatrix(int[],int,int)... ",
                  "submatrix not successfully retreived")
  except IndexDefect:
    tryFailure(errorCount, "getMatrix(int[],int,int)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  try:
    M = B[badrowindexset, columnindexset]
    tryFailure(errorCount, "getMatrix(int[],int[])... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      M = B[rowindexset, badcolumnindexset]
      tryFailure(errorCount, "getMatrix(int[],int[])... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("getMatrix(int[],int[])... ArrayIndexOutOfBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "getMatrix(int[],int[])... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    M = B[rowindexset, columnindexset]
    try:
      check(SUB, M)
      trySuccess("getMatrix(int[],int[])... ", "")
    except ValueError:
      tryFailure(errorCount, "getMatrix(int[],int[])... ",
                  "submatrix not successfully retreived")
  except IndexDefect:
    tryFailure(errorCount, "getMatrix(int[],int[])... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  # Various set methods:
  try:
    B[B.rowDimension, B.columnDimension-1] = 0
    tryFailure(errorCount, "set(int,int,double)... ",
                "OutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      B[B.rowDimension-1, B.columnDimension] = 0
      tryFailure(errorCount, "set(int,int,double)... ",
                  "OutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("set(int,int,double)... OutofBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "set(int,int,double)... ",
                "OutOfBoundsException expected but not thrown")
  try:
    B[ib, jb] = 0
    tmp = B[ib, jb]
    try:
      check(tmp, 0)
      trySuccess("set(int,int,double)... ", "")
    except ValueError:
      tryFailure(errorCount, "set(int,int,double)... ",
                  "Matrix element not successfully set")
  except IndexDefect:
    tryFailure(errorCount, "set(int,int,double)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  M = matrix[float32](2, 3)
  try:
    B[ib .. ie+B.rowDimension+1, jb .. je] = M
    tryFailure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      B[ib .. ie, jb .. je+B.columnDimension+1] = M
      tryFailure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("setMatrix(int,int,int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "")
  except AssertionDefect:
    tryFailure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    B[ib .. ie, jb .. je] = M
    try:
      check(M - B[ib .. ie, jb .. je], M)
      trySuccess("setMatrix(int,int,int,int,Matrix)... ", "")
    except ValueError:
      tryFailure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                  "submatrix not successfully set")
    B[ib .. ie, jb .. je] = SUB
  except IndexDefect:
    tryFailure(errorCount, "setMatrix(int,int,int,int,Matrix)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  try:
    B[ib .. ie+B.rowDimension+1, columnindexset] = M
    tryFailure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      B[ib .. ie, badcolumnindexset] = M
      tryFailure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("setMatrix(int,int,int[],Matrix)... ArrayIndexOutOfBoundsException... ", "")

  except AssertionDefect:
    tryFailure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")

  try:
    B[ib .. ie, columnindexset] = M
    try:
      check(M - B[ib .. ie, columnindexset], M)
      trySuccess("setMatrix(int,int,int[],Matrix)... ", "")
    except ValueError:
      tryFailure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                  "submatrix not successfully set")
    B[ib .. ie, jb .. je] = SUB
  except IndexDefect:
    tryFailure(errorCount, "setMatrix(int,int,int[],Matrix)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  try:
    B[rowindexset, jb .. je+B.columnDimension+1] = M
    tryFailure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      B[badrowindexset, jb .. je] = M
      tryFailure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("setMatrix(int[],int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "")

  except AssertionDefect:
    tryFailure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    B[rowindexset, jb .. je] = M
    try:
      check(M - B[rowindexset, jb .. je], M)
      trySuccess("setMatrix(int[],int,int,Matrix)... ", "")
    except ValueError:
      tryFailure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                  "submatrix not successfully set")
    B[ib .. ie, jb .. je] = SUB
  except IndexDefect:
    tryFailure(errorCount, "setMatrix(int[],int,int,Matrix)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
  try:
    B[rowindexset, badcolumnindexset] = M
    tryFailure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  except IndexDefect:
    try:
      B[badrowindexset, columnindexset] = M
      tryFailure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
    except IndexDefect:
      trySuccess("setMatrix(int[],int[],Matrix)... ArrayIndexOutOfBoundsException... ", "")

  except AssertionDefect:
    tryFailure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                "ArrayIndexOutOfBoundsException expected but not thrown")
  try:
    B[rowindexset, columnindexset] = M
    try:
      check(M - B[rowindexset, columnindexset], M)
      trySuccess("setMatrix(int[],int[],Matrix)... ", "")
    except ValueError:
      tryFailure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                  "submatrix not successfully set")
  except IndexDefect:
    tryFailure(errorCount, "setMatrix(int[],int[],Matrix)... ",
                "Unexpected ArrayIndexOutOfBoundsException")
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
  A = matrix(columnwise, 3)
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

  echo("\nTestMatrix completed.")
  echo("Total errors reported: ", errorCount)
  #echo("Total warnings reported: ", warningCount)
  if errorCount > 0: raise newException(ValueError, "")

try:
  main()
except:
  quit(QuitFailure)
