# Morpheus — Nim Matrix library

This is a port of the NIST JAMA library to Nim.

API documentation is [here](https://rawgit.com/notTito/morpheus/master/htmldocs/morpheus.html)

## Example of use

```nim
   import morpheus
   # Solve a linear system A x = b and compute the residual norm, ||b - A x||.
   let vals = @[@[1.0, 2.0, 3.0], @[4.0, 5.0, 6.0], @[7.0, 8.0, 10.0]]
   let A = matrix(vals)
   let b = randMatrix(3, 1)
   let x = A.solve(b)
   let r = A * x - b
   let rnorm = r.normInf()
```

## Todo
- Port the EigenValueDecomposition class

## Copyright Notice
From the original JAMA code:

> This software is a cooperative product of The MathWorks and the National 
> Institute of Standards and Technology (NIST) which has been released to the 
> public domain. Neither The MathWorks nor NIST assumes any responsibility 
> whatsoever for its use by other parties, and makes no guarantees, expressed or 
> implied, about its quality, reliability, or any other characteristic.
