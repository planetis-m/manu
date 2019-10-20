# See: https://peterroelants.github.io/posts/neural-network-implementation-part04/
# Implementing an Artificial Neural Network in manu
import math, strutils, "../manu" / matrix

proc sigmoid(s: float): float {.inline.} =
   result = 1.0 / (1.0 + exp(-s))
makeUniversal(sigmoid)

proc loss(y, t: float): float {.inline.} =
   result = t * ln(y) + (1.0 - t) * ln(1.0 - y)
makeUniversalBinary(loss)

proc main =
   const
      m = 4 # batch length
      nodes = 3
      rate = 0.01
   let
      X = matrix(2, @[0.0, 0, 0, 1, 1, 0, 1, 1])
      Y = matrix(1, @[0.0, 1, 1, 0])
   var
      # LAYER 1
      W1 = randMatrix(2, nodes, -1.0..1.0)
      b1 = rowVector(nodes)
      # LAYER 2
      W2 = randMatrix(nodes, 1, -1.0..1.0)
      b2 = 0.0
   for i in 1 .. 1000:
      # Foward Prop
      let
         # LAYER 1
         Z1 = X * W1 + b1 # broadcast bias to (m, nodes)
         A1 = sigmoid(Z1)
         # LAYER 2
         Z2 = A1 * W2 + b2 # to (m, 1)
         A2 = sigmoid(Z2)
      # Back Prop
      let
         # LAYER 2
         dZ2 = A2 - Y
         db2 = sum(dZ2)
         dW2 = A1.transpose * dZ2
         # LAYER 1
         dZ1 = (dZ2 * W2.transpose) *. (1.0 - A1) *. A1
         db1 = sumColumns(dZ1)
         dW1 = X.transpose * dZ1
      # Gradient Descent
      # LAYER 2
      W2 -= rate * dW2
      b2 -= rate * db2
      # LAYER 1
      W1 -= rate * dW1
      b1 -= rate * db1
      # Cross Entropy
      let loss = -sum(loss(A2, Y)) / m.float
      if i mod 20 == 0:
         echo(" Iteration ", i, ":")
         echo("   Loss = ", formatEng(loss))
         echo("   Predictions =\n", A2)

main()
