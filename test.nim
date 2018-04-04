import morpheus

const vals = @[@[1.0, 2.0, 3.0], @[4.0, 5.0, 6.0], @[7.0, 8.0, 10.0]]
let a = matrix(vals)
let b = randMatrix(3, 1)
let x = a.solve(b)
let r = a * x - b
let rnorm = r.normInf()

echo x
