import matrix, cholesky, qr, lu#, eigen, svd
export matrix, cholesky, qr, lu

# proc norm2*(m: Matrix): float =
#    ## Two norm,
#    ## returns maximum singular value.
#    svd(m).norm2())

proc solve*(a, b: Matrix): Matrix =
   ## Solve ``A*X = B``,
   ## returns solution if A is square, least squares solution otherwise
   ## parameter ``b``: the right hand side
   if a.m == a.n:
      lu(a).solve(b)
   else:
      qr(a).solve(b)

proc solveTranspose*(a, b: Matrix): Matrix =
   ## Solve ``X*A = B``, which is also ``A'*X' = B'``,
   ## returns solution if A is square, least squares solution otherwise.
   ## parameter ``b``: the right hand side
   transpose(a).solve(b.transpose())

proc inverse*(m: Matrix): Matrix =
   ## Matrix inverse or pseudoinverse,
   ## returns inverse(A) if A is square, pseudoinverse otherwise.
   solve(m, identity(m.m, m.m))

proc det*(m: Matrix): float =
   ## Matrix determinant
   lu(m).det()

# proc rank*(m: Matrix): int =
#    ## Matrix rank,
#    ## returns effective numerical rank, obtained from SVD.
#    svd(m).rank()

# proc cond*(m: Matrix): float =
#    ## Matrix condition (2 norm),
#    ## returns ratio of largest to smallest singular value.
#    svd(m).cond()
