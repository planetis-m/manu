import
   math,
   testutils,
   morpheus / [matrix, cholesky]

block:
   let
      columnwise = @[1.0, -2, -2, 5]
      mat = matrix(columnwise, 2)
      l = chol(mat).getL

   check(l * l.transpose, mat)

block:
   let
      columnwise = @[25.0, 15, -5, 15, 18, 0, -5, 0, 11]
      mat = matrix(columnwise, 3)
      l = chol(mat).getL

   check(l * l.transpose, mat)

block:
   let
      columnwise = @[18.0, 22, 54, 42, 22, 70, 86, 62, 54, 86, 174, 134, 42, 62, 134, 106]
      mat = matrix(columnwise, 4)
      l = chol(mat).getL

   check(l * l.transpose, mat)

block:
   let
      columnwise = @[6.0, 15, 55, 72, 101,
                     15, 55, 225, 229, 256,
                     55, 225, 979, 1024, 1200,
                     72, 229, 1024, 2048, 2057,
                     101, 256, 1200, 2057, 6000]
      mat = matrix(columnwise, 5)
      l = chol(mat).getL

   check(l * l.transpose, mat)
