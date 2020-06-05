import random, strutils, "../manu" / matrix

proc loss(y, t: float): float {.inline.} =
   0.5 * (y - t) * (y - t)
makeUniversalBinary(loss)

proc celsiusToFahrenheit(c: float): float {.inline.} =
   result = c * 1.8 + 32
makeUniversal(celsiusToFahrenheit)

proc main =
   const
      m = 100
      epochs = 7_000
      rate = 0.0001
   let
      X = randMatrix(1, m, 1.0..100.0)
      Y = celsiusToFahrenheit(X)
   var
      # LAYER 1
      W = rand(-1.0..1.0)
      b = 0.0
   for i in 1 .. epochs:
      # Foward Prop
      let
         Z = X * W + b
      # Back Prop
      let
         dZ = Z - Y
         db = sum(dZ) / m.float
         dW = sum(X.transpose * dZ) / m.float
      # Gradient Descent
      W -= rate * dW
      b -= rate * db
      let loss = sum(loss(Z, Y)) / m.float
      if i mod 250 == 0:
         echo(" Iteration ", i, ":")
         echo("   Loss = ", formatEng(loss))

main()
