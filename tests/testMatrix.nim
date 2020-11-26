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
import math, "../manu", testutils

# main

proc main() =
   var A, B, C, Z, O, I, R, S, X, SUB, M, T, SQ, DEF, SOL: Matrix[float32]
   var errorCount = 0
   var warningCount = 0
   var tmp, s: float
   let columnwise = @[1.0'f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
   let rowwise = @[1.0'f32, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 3.0, 6.0, 9.0, 12.0]
   let avals = @[@[1.0'f32, 4.0, 7.0, 10.0], @[2.0'f32, 5.0, 8.0, 11.0], @[3.0'f32, 6.0, 9.0, 12.0]]
   let rankdef = avals
   let tvals =  @[@[1.0'f32, 2.0, 3.0], @[4.0'f32, 5.0, 6.0], @[7.0'f32, 8.0, 9.0], @[10.0'f32, 11.0, 12.0]]
   let subavals = @[@[5.0'f32, 8.0, 11.0], @[6.0'f32, 9.0, 12.0]]
   var rvals = @[@[2.0'f32, 5.0, 8.0, 11.0], @[3.0'f32, 6.0, 9.0, 12.0]]
   rvals.add @[1.0'f32, 4.0, 7.0]
   let pvals = @[@[4.0'f32, 1.0, 1.0], @[1.0'f32, 2.0, 3.0], @[1.0'f32, 3.0, 6.0]]
   let ivals = @[@[1.0'f32, 0.0, 0.0], @[0.0'f32, 1.0, 0.0], @[0.0'f32, 0.0, 1.0]]
   let evals =
      @[@[0.0'f32, 1.0, 0.0, 0.0], @[1.0'f32, 0.0, 2.0e-7,0.0], @[0.0'f32, -2.0e-7,0.0, 1.0], @[0.0'f32, 0.0, 1.0, 0.0]]
   let square = @[@[166.0'f32, 188.0, 210.0], @[188.0'f32, 214.0, 240.0], @[210.0'f32, 240.0, 270.0]]
   let sqSolution = @[@[13.0'f32], @[15.0'f32]]
   let condmat = @[@[1.0'f32, 3.0], @[7.0'f32, 9.0]]
   let badeigs = @[@[0.0'f32, 0, 0, 0, 0], @[0.0'f32, 0, 0, 0, 1], @[0.0'f32, 0, 0, 1, 0],
      @[1.0'f32, 1, 0, 0, 1], @[1.0'f32, 0, 1, 0, 1]]
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
   let columnsummax = 33.0
   let rowsummax = 30.0
   let sumofdiagonals = 15.0
   let sumofsquares = 650.0

   # Constructors
   echo("Testing constructors and constructor-like methods...")
   try:
      # check that exception is thrown in packed constructor with invalid length
      A = matrix(columnwise, invalidld)
      try_failure(errorCount, "Catch invalid length in packed constructor... ",
                  "exception not thrown for invalid input")
   except AssertionDefect as e:
      try_success("Catch invalid length in packed constructor... ", e.msg)
   try:
      # check that exception is thrown in default constructor if input array is 'ragged'
      A = matrix(rvals)
      tmp = A[raggedr, raggedc]
   except AssertionDefect as e:
      try_success("Catch ragged input to default constructor... ", e.msg)
   except IndexDefect:
      try_failure(errorCount, "Catch ragged input to constructor... ",
                  "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later")
#    try:
#       # check that exception is thrown in constructWithCopy if input array is 'ragged'
#       A = constructWithCopy(rvals)
#       tmp = A[raggedr, raggedc]
#    except AssertionDefect as e:
#       try_success("Catch ragged input to constructWithCopy... ", e.msg)
#    except IndexDefect:
#       try_failure(errorCount, "Catch ragged input to constructWithCopy... ",
#                   "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later")

   A = matrix(columnwise, validld)
   B = matrix(avals)
#    tmp = B[0, 0]
#    avals[0][0] = 0.0
   C = B - A
#    avals[0][0] = tmp
#    B = constructWithCopy(avals)
#    tmp = B[0, 0]
#    avals[0][0] = 0.0
#    if tmp - B[0,0] != 0.0:
#       # check that constructWithCopy behaves properly
#       try_failure(errorCount, "constructWithCopy... ",
#                   "copy not effected... data visible outside")
#    else:
#       try_success("constructWithCopy... ", "")
#
#    avals[0][0] = columnwise[0]
   I = matrix(ivals)
   try:
      check(I, identity[float32](3))
      try_success("identity... ", "")
   except ValueError:
      try_failure(errorCount, "identity... ",
                  "identity Matrix not successfully created")

   # Access Methods:
   echo("\nTesting access methods...")

   # Various get methods:
   B = matrix(avals)
   if B.rowDimension != rows:
      try_failure(errorCount,"rowDimension... ", "")
   else:
      try_success("getRowDimension... ", "")

   if B.columnDimension != cols:
      try_failure(errorCount,"columnDimension... ", "")
   else:
      try_success("ColumnDimension... ", "")
   B = matrix(avals)
   var barray = B.getArray()
   if barray != avals:
      try_failure(errorCount,"getArray... ", "")
   else:
      try_success("getArray... ", "")
#    barray = B.getArrayCopy()
#    if barray == avals:
#       try_failure(errorCount,"getArrayCopy... ", "data not (deep) copied")
#    try:
#       check(barray, avals)
#       try_success("getArrayCopy... ", "")
#    except ValueError:
#       try_failure(errorCount,"getArrayCopy... ","data not successfully (deep) copied")

   var bpacked = B.getColumnPacked()
   try:
      check(bpacked, columnwise)
      try_success("getColumnPacked... ", "")
   except ValueError:
      try_failure(errorCount,"getColumnPacked... ",
                  "data not successfully (deep) copied by columns")
   bpacked = B.getRowPacked()
   try:
      check(bpacked, rowwise)
      try_success("getRowPacked... ", "")
   except ValueError:
      try_failure(errorCount,"getRowPacked... ",
                  "data not successfully (deep) copied by rows")
   try:
      tmp = B[B.rowDimension, B.columnDimension-1]
      try_failure(errorCount, "get(int,int)... ",
                  "OutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         tmp = B[B.rowDimension-1, B.columnDimension]
         try_failure(errorCount,"get(int,int)... ",
                     "OutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("get(int,int)... OutofBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"get(int,int)... ",
                  "OutOfBoundsException expected but not thrown")
   try:
      if B[B.rowDimension-1, B.columnDimension-1] !=
            avals[B.rowDimension-1][B.columnDimension-1]:
         try_failure(errorCount,"get(int,int)... ",
                     "Matrix entry (i,j) not successfully retreived")
      else:
         try_success("get(int,int)... ", "")
   except IndexDefect:
      try_failure(errorCount,"get(int,int)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   SUB = matrix(subavals)
   try:
      M = B[ib .. ie+B.rowDimension+1, jb .. je]
      try_failure(errorCount,"getMatrix(int,int,int,int)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         M = B[ib .. ie, jb .. je+B.columnDimension+1]
         try_failure(errorCount,"getMatrix(int,int,int,int)... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("getMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"getMatrix(int,int,int,int)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      M = B[ib .. ie, jb .. je]
      try:
         check(SUB, M)
         try_success("getMatrix(int,int,int,int)... ", "")
      except ValueError:
         try_failure(errorCount,"getMatrix(int,int,int,int)... ",
                     "submatrix not successfully retreived")
   except IndexDefect:
      try_failure(errorCount,"getMatrix(int,int,int,int)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   try:
      M = B[ib .. ie, badcolumnindexset]
      try_failure(errorCount,"getMatrix(int,int,int[])... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         M = B[ib .. ie+B.rowDimension+1, columnindexset]
         try_failure(errorCount,"getMatrix(int,int,int[])... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("getMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"getMatrix(int,int,int[])... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      M = B[ib .. ie, columnindexset]
      try:
         check(SUB, M)
         try_success("getMatrix(int,int,int[])... ", "")
      except ValueError:
         try_failure(errorCount,"getMatrix(int,int,int[])... ",
                     "submatrix not successfully retreived")
   except IndexDefect:
      try_failure(errorCount,"getMatrix(int,int,int[])... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   try:
      M = B[badrowindexset, jb .. je]
      try_failure(errorCount,"getMatrix(int[],int,int)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         M = B[rowindexset, jb .. je+B.columnDimension+1]
         try_failure(errorCount,"getMatrix(int[],int,int)... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("getMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"getMatrix(int[],int,int)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      M = B[rowindexset, jb .. je]
      try:
         check(SUB, M)
         try_success("getMatrix(int[],int,int)... ", "")
      except ValueError:
         try_failure(errorCount,"getMatrix(int[],int,int)... ",
                     "submatrix not successfully retreived")
   except IndexDefect:
      try_failure(errorCount,"getMatrix(int[],int,int)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   try:
      M = B[badrowindexset, columnindexset]
      try_failure(errorCount,"getMatrix(int[],int[])... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         M = B[rowindexset, badcolumnindexset]
         try_failure(errorCount,"getMatrix(int[],int[])... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("getMatrix(int[],int[])... ArrayIndexOutOfBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"getMatrix(int[],int[])... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      M = B[rowindexset, columnindexset]
      try:
         check(SUB, M)
         try_success("getMatrix(int[],int[])... ", "")
      except ValueError:
         try_failure(errorCount,"getMatrix(int[],int[])... ",
                     "submatrix not successfully retreived")
   except IndexDefect:
      try_failure(errorCount,"getMatrix(int[],int[])... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   # Various set methods:
   try:
      B[B.rowDimension,B.columnDimension-1] = 0.0
      try_failure(errorCount,"set(int,int,double)... ",
                  "OutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         B[B.rowDimension-1,B.columnDimension] = 0.0
         try_failure(errorCount,"set(int,int,double)... ",
                     "OutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("set(int,int,double)... OutofBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"set(int,int,double)... ",
                  "OutOfBoundsException expected but not thrown")
   try:
      B[ib, jb] = 0.0
      tmp = B[ib, jb]
      try:
         check(tmp, 0.0)
         try_success("set(int,int,double)... ", "")
      except ValueError:
         try_failure(errorCount,"set(int,int,double)... ",
                     "Matrix element not successfully set")
   except IndexDefect:
      try_failure(errorCount,"set(int,int,double)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   M = matrix[float32](2, 3)
   try:
      B[ib .. ie+B.rowDimension+1, jb .. je] = M
      try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         B[ib .. ie, jb .. je+B.columnDimension+1] = M
         try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("setMatrix(int,int,int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "")
   except AssertionDefect:
      try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      B[ib .. ie, jb .. je] = M
      try:
         check(M - B[ib .. ie, jb .. je], M)
         try_success("setMatrix(int,int,int,int,Matrix)... ", "")
      except ValueError:
         try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ",
                     "submatrix not successfully set")
      B[ib .. ie, jb .. je] = SUB
   except IndexDefect:
      try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   try:
      B[ib .. ie+B.rowDimension+1, columnindexset] = M
      try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         B[ib .. ie, badcolumnindexset] = M
         try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("setMatrix(int,int,int[],Matrix)... ArrayIndexOutOfBoundsException... ", "")

   except AssertionDefect:
      try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")

   try:
      B[ib .. ie, columnindexset] = M
      try:
         check(M - B[ib .. ie, columnindexset], M)
         try_success("setMatrix(int,int,int[],Matrix)... ", "")
      except ValueError:
         try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ",
                     "submatrix not successfully set")
      B[ib .. ie, jb .. je] = SUB
   except IndexDefect:
      try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   try:
      B[rowindexset, jb .. je+B.columnDimension+1] = M
      try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         B[badrowindexset, jb .. je] = M
         try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("setMatrix(int[],int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "")

   except AssertionDefect:
      try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      B[rowindexset, jb .. je] = M
      try:
         check(M - B[rowindexset, jb .. je], M)
         try_success("setMatrix(int[],int,int,Matrix)... ", "")
      except ValueError:
         try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ",
                     "submatrix not successfully set")
      B[ib .. ie, jb .. je] = SUB
   except IndexDefect:
      try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   try:
      B[rowindexset, badcolumnindexset] = M
      try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   except IndexDefect:
      try:
         B[badrowindexset, columnindexset] = M
         try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ",
                     "ArrayIndexOutOfBoundsException expected but not thrown")
      except IndexDefect:
         try_success("setMatrix(int[],int[],Matrix)... ArrayIndexOutOfBoundsException... ", "")

   except AssertionDefect:
      try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ",
                  "ArrayIndexOutOfBoundsException expected but not thrown")
   try:
      B[rowindexset, columnindexset] = M
      try:
         check(M - B[rowindexset, columnindexset], M)
         try_success("setMatrix(int[],int[],Matrix)... ", "")
      except ValueError:
         try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ",
                     "submatrix not successfully set")
   except IndexDefect:
      try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ",
                  "Unexpected ArrayIndexOutOfBoundsException")
   # Array-like methods:
   echo("\nTesting array-like methods...")
   S = matrix(columnwise,nonconformld)
   R = randMatrix[float32](A.rowDimension, A.columnDimension)
   A = R
   try:
      S = A - S
      try_failure(errorCount,"minus conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("minus conformance check... ", "")
   if (A - R).norm1() != 0.0:
      try_failure(errorCount,"minus... ",
                  "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)")
   else:
      try_success("minus... ", "")
   A = R #.copy()
   A -= R
   Z = matrix[float32](A.rowDimension,A.columnDimension)
   try:
      A -= S
      try_failure(errorCount,"minusEquals conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("minusEquals conformance check... ", "")
   if (A - Z).norm1() != 0.0:
      try_failure(errorCount,"minusEquals... ",
                  "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)")
   else:
      try_success("minusEquals... ", "")
   A = R #.copy()
   B = randMatrix[float32](A.rowDimension, A.columnDimension)
   C = A - B
   try:
      S = A + S
      try_failure(errorCount,"plus conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("plus conformance check... ", "")
   try:
      check(C + B, A)
      try_success("plus... ", "")
   except ValueError:
      try_failure(errorCount,"plus... ","(C = A - B, but C + B != A)")
   C = A - B
   C += B
   try:
      A += S
      try_failure(errorCount,"plusEquals conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("plusEquals conformance check... ", "")
   try:
      check(C, A)
      try_success("plusEquals... ", "")
   except ValueError:
      try_failure(errorCount,"plusEquals... ","(C = A - B, but C = C + B != A)")
   A = -R
   try:
      check(A + R, Z)
      try_success("uminus... ", "")
   except ValueError:
      try_failure(errorCount,"uminus... ","(-A + A != zeros)")
   A = R #.copy()
   O = matrix(A.rowDimension,A.columnDimension, 1.0'f32)
   #C = A \. R
   #try:
      #S = A \. S
      #try_failure(errorCount,"arrayLeftDivide conformance check... ","nonconformance not raised")
   #except AssertionDefect:
      #try_success("arrayLeftDivide conformance check... ", "")
   #try:
      #check(C, O)
      #try_success("arrayLeftDivide... ", "")
   #except ValueError:
      #try_failure(errorCount,"arrayLeftDivide... ","(M.\\M != ones)")
   #try:
      #A \.= S
      #try_failure(errorCount,"arrayLeftDivideEquals conformance check... ","nonconformance not raised")
   #except AssertionDefect:
      #try_success("arrayLeftDivideEquals conformance check... ", "")
   #A \.= R
   #try:
      #check(A,O)
      #try_success("arrayLeftDivideEquals... ", "")
   #except ValueError:
      #try_failure(errorCount,"arrayLeftDivideEquals... ","(M.\\M != ones)")
   #A = R #.copy()
   try:
      C = A /. S
      try_failure(errorCount,"arrayRightDivide conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("arrayRightDivide conformance check... ", "")
   C = A /. R
   try:
      check(C, O)
      try_success("arrayRightDivide... ", "")
   except ValueError:
      try_failure(errorCount,"arrayRightDivide... ","(M./M != ones)")
   try:
      A /.= S
      try_failure(errorCount,"arrayRightDivideEquals conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("arrayRightDivideEquals conformance check... ", "")
   A /.= R
   try:
      check(A, O)
      try_success("arrayRightDivideEquals... ", "")
   except ValueError:
      try_failure(errorCount,"arrayRightDivideEquals... ","(M./M != ones)")
   A = R #.copy()
   B = randMatrix[float32](A.rowDimension, A.columnDimension)
   try:
      S = A *. S
      try_failure(errorCount, "arrayTimes conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("arrayTimes conformance check... ", "")
   C = A *. B
   try:
      check(C /. B, A)
      try_success("arrayTimes... ", "")
   except ValueError:
      try_failure(errorCount,"arrayTimes... ","(A = R, C = A.*B, but C./B != A)")
   try:
      A *.= S
      try_failure(errorCount,"arrayTimesEquals conformance check... ","nonconformance not raised")
   except AssertionDefect:
      try_success("arrayTimesEquals conformance check... ", "")
   A *.= B
   try:
      check(A /. B, R)
      try_success("arrayTimesEquals... ", "")
   except ValueError:
      try_failure(errorCount,"arrayTimesEquals... ","(A = R, A = A.*B, but A./B != R)")

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
      try_success("transpose...", "")
   except ValueError:
      try_failure(errorCount,"transpose()...","transpose unsuccessful")
   try:
      check(A.norm1(), columnsummax)
      try_success("norm1...", "")
   except ValueError:
      try_failure(errorCount,"norm1()...", "incorrect norm calculation")
   try:
      check(A.normInf(),rowsummax)
      try_success("normInf()...", "")
   except ValueError:
      try_failure(errorCount,"normInf()...","incorrect norm calculation")
   try:
      check(A.normF(), sqrt(sumofsquares))
      try_success("normF...", "")
   except ValueError:
      try_failure(errorCount,"normF()...","incorrect norm calculation")
   try:
      check(A.trace(), sumofdiagonals)
      try_success("trace()...", "")
   except ValueError:
      try_failure(errorCount,"trace()...","incorrect trace calculation")
   try:
      check(A[0 .. A.rowDimension-1, 0 .. A.rowDimension-1].det(), 0.0)
      try_success("det()...", "")
   except ValueError:
      try_failure(errorCount,"det()...",
                  "incorrect determinant calculation")
   SQ = matrix(square)
   try:
      check(A * A.transpose(), SQ)
      try_success("times(Matrix)...", "")
   except ValueError:
      try_failure(errorCount,"times(Matrix)...",
                  "incorrect Matrix-Matrix product calculation")
   try:
      check(A * 0.0, Z)
      try_success("times(double)...", "")
   except ValueError:
      try_failure(errorCount,"times(double)...",
                  "incorrect Matrix-scalar product calculation")
   A = matrix(columnwise, 4)
   let QR = A.qr()
   R = QR.getR()
   try:
      check(A, QR.getQ() * R)
      try_success("QRDecomposition...", "")
   except ValueError:
      try_failure(errorCount,"QRDecomposition...","incorrect QR decomposition calculation")
   var SVD = A.svd()
   try:
      check(A, SVD.getU() * SVD.getS() * SVD.getV().transpose())
      try_success("SingularValueDecomposition...", "")
   except ValueError:
      try_failure(errorCount,"SingularValueDecomposition...",
                  "incorrect singular value decomposition calculation")
   DEF = matrix(rankdef)
   try:
      assert(DEF.rank() == min(DEF.rowDimension,DEF.columnDimension)-1)
      try_success("rank()...", "")
   except AssertionDefect:
      try_failure(errorCount,"rank()...","incorrect rank calculation")
   B = matrix(condmat)
   SVD = B.svd()
   let singularvalues = SVD.getSingularValues()
   try:
      check(B.cond(), singularvalues[0] / singularvalues[min(B.rowDimension,B.columnDimension)-1])
      try_success("cond()...", "")
   except ValueError:
      try_failure(errorCount,"cond()...","incorrect condition number calculation")
   let n = A.columnDimension
   A = A[0 .. n-1, 0 .. n-1]
   A[0, 0] = 0.0
   let LU = A.lu()
   try:
      check(A[LU.getPivot(), 0 .. n-1], LU.getL() * LU.getU())
      try_success("LUDecomposition...", "")
   except ValueError:
      try_failure(errorCount,"LUDecomposition...","incorrect LU decomposition calculation")
   X = A.inverse()
   try:
      check(A * X, identity[float32](3))
      try_success("inverse()...", "")
   except ValueError:
      try_failure(errorCount,"inverse()...","incorrect inverse calculation")
   O = matrix(SUB.rowDimension,1,1.0'f32)
   SOL = matrix(sqSolution)
   SQ = SUB[0 .. SUB.rowDimension-1, 0 .. SUB.rowDimension-1]
   try:
      check(SQ.solve(SOL), O)
      try_success("solve()...", "")
   except AssertionDefect as e1:
      try_failure(errorCount,"solve()...", e1.msg)
   except ValueError as e:
      try_failure(errorCount,"solve()...", e.msg)
   A = matrix(pvals)
   let Chol = A.chol()
   let L = Chol.getL()
   try:
      check(A, L * L.transpose())
      try_success("CholeskyDecomposition...", "")
   except ValueError:
      try_failure(errorCount,"CholeskyDecomposition...",
                  "incorrect Cholesky decomposition calculation")
   X = Chol.solve(identity[float32](3))
   try:
      check(A * X, identity[float32](3))
      try_success("CholeskyDecomposition solve()...", "")
   except ValueError:
      try_failure(errorCount,"CholeskyDecomposition solve()...",
                  "incorrect Choleskydecomposition solve calculation")
   var Eig = A.eig()
   var D = Eig.getD()
   var V = Eig.getV()
   try:
      check(A * V, V * D)
      try_success("EigenvalueDecomposition (symmetric)...", "")
   except ValueError:
      try_failure(errorCount, "EigenvalueDecomposition (symmetric)...",
                  "incorrect symmetric Eigenvalue decomposition calculation")
   A = matrix(evals)
   Eig = A.eig()
   D = Eig.getD()
   V = Eig.getV()
   try:
      check(A * V, V * D)
      try_success("EigenvalueDecomposition (nonsymmetric)...", "")
   except ValueError:
      try_failure(errorCount,"EigenvalueDecomposition (nonsymmetric)...",
                  "incorrect nonsymmetric Eigenvalue decomposition calculation")
   try:
      echo("\nTesting Eigenvalue; If this hangs, we've failed")
      let bA = matrix(badeigs)
      let bEig = bA.eig()
      try_success("EigenvalueDecomposition (hang)...", "")
   except ValueError:
      try_failure(errorCount,"EigenvalueDecomposition (hang)...",
                  "incorrect termination")

   echo("\nTestMatrix completed.")
   echo("Total errors reported: ", errorCount)
   #echo("Total warnings reported: ", warningCount)
   if errorCount > 0: quit(QuitFailure)

main()
