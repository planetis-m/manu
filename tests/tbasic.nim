import manu, utils

# main
proc main() =
  var A, B, C, I, SUB, M: Matrix[float32]
  var errorCount = 0
  var warningCount = 0
  var tmp: float32
  let columnwise = @[1'f32, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
  let rowwise = @[1'f32, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12]
  let avals = @[@[1'f32, 4, 7, 10], @[2'f32, 5, 8, 11], @[3'f32, 6, 9, 12]]
  let subavals = @[@[5'f32, 8, 11], @[6'f32, 9, 12]]
  var rvals = @[@[2'f32, 5, 8, 11], @[3'f32, 6, 9, 12]]
  rvals.add @[1'f32, 4, 7]
  let ivals = @[@[1'f32, 0, 0], @[0'f32, 1, 0], @[0'f32, 0, 1]]
  let
    rows = 3
    cols = 4
  let invalidld = 5 # should trigger bad shape for construction with val
  let
    raggedr = 0 # (raggedr,raggedc) should be out of bounds in ragged array
    raggedc = 4
  let validld = 3 # leading dimension of intended test Matrices
  let
    ib = 1
    ie = 2
    jb = 1
    je = 3 # index ranges for sub Matrix
  let rowindexset = [1, 2]
  let badrowindexset = [1, 3]
  let columnindexset = [1, 2, 3]
  let badcolumnindexset = [1, 2, 4]

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

  echo("\nTestBasic completed.")
  echo("Total errors reported: ", errorCount)
  #echo("Total warnings reported: ", warningCount)
  if errorCount > 0: raise newException(ValueError, "")

try:
  main()
except:
  quit(QuitFailure)
