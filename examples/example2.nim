import manu

# Solve a linear system A x = b and compute the residual norm, ||b - A x||.
let vals = @[@[1.0, 2.0, 3.0], @[4.0, 5.0, 6.0], @[7.0, 8.0, 10.0]]
let A = matrix(vals)
let b = randMatrix(3, 1)
let x = A.solve(b)
let r = A * x - b
let rnorm = r.normInf()
echo("x =\n", x)
echo("residual norm = ", rnorm)
