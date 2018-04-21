import
   math,
   testutils,
   morpheus / [matrix, lu]

block:
   let
      columnwise = @[14.0, 6, 7, 18]
      mat = matrix(columnwise, 2)
      lu = mat.lu()

   check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)

block:
   let
      columnwise = @[1.0, 0, 2, 0, 10, 0, 2, 0, 9]
      mat = matrix(columnwise, 3)
      lu = mat.lu()

   check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)

block:
   let
      columnwise = @[77.0, -9.0, 11.0, 19.0, 17.0, 100.0, 1.0, 34.0, -2.0]
      mat = matrix(columnwise, 3)
      lu = mat.lu()

   check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)

block:
   let
      columnwise = @[99.0, 1, -10, 6, 14, 65, 7, 48, 39, 40, -2, 9, 11, 5, 43, 99]
      mat = matrix(columnwise, 4)
      lu = mat.lu()

   check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
