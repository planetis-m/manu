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
      columnwise = @[77.0, 19, 1, -9, 17, 34, 11, 100, -2]
      l = matrix(@[1.0, 0, 0, 0.143, 1, 0, -0.117, 0.198, 1], 3)
      u = matrix(@[77.0, 19, 1, 0, 97.286, -2.143, 0, 0, 34.540], 3)
      piv = @[1, 0, 0, 0, 0, 1, 0, 1, 0]
      mat = matrix(columnwise, 3)
      lu = mat.lu()

   #check(lu.getL, l)
   #check(lu.getU, u)
   #assert lu.getPivot == piv
   check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)

block:
   let
      columnwise = @[99.0, 1, -10, 6, 14, 65, 7, 48, 39, 40, -2, 9, 11, 5, 43, 99]
      mat = matrix(columnwise, 4)
      lu = mat.lu()

   check(mat[lu.getPivot, 0 ..< mat.n], lu.getL * lu.getU)
