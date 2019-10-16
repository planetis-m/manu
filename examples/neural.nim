# From: https://gist.github.com/Jeraldy/1aa6ae6fefa46b7a9cc02b6573cfeefe
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
      nodes = 400
      anneal = 0.99
   let
      X = matrix(m, @[0.0, 0, 1, 1, 0, 1, 0, 1]) # X.transpose
      Y = matrix(m, @[0.0, 1, 1, 0]) # Y.transpose
   var
      rate = 1.0
      # LAYER 1
      W1 = randMatrix(nodes, 2)
      b1 = zeros(nodes, m)
      # LAYER 2
      W2 = randMatrix(1, nodes)
      b2 = zeros(1, m)
   for i in 1 .. 100:
      # Foward Prop
      let
         # LAYER 1
         Z1 = W1 * X + b1
         A1 = sigmoid(Z1)
         # LAYER 2
         Z2 = W2 * A1 + b2
         A2 = sigmoid(Z2)
      # Back Prop
      let
         # LAYER 2
         dZ2 = A2 - Y
         db2 = dZ2 / m.float
         dW2 = db2 * A1.transpose
         # LAYER 1
         dZ1 = (W2.transpose * dZ2) *. (1.0 - A1 *. A1)
         db1 = dZ1 / m.float
         dW1 = db1 * X.transpose
      # Gradient Descent
      # LAYER 2
      W2 -= rate * dW2
      b2 -= rate * db2
      # LAYER 1
      W1 -= rate * dW1
      b1 -= rate * db1
      rate = rate * anneal
      # Cross Entropy
      let loss = -sum(loss(A2, Y)) / m.float
      if i mod 20 == 0:
         echo(" Iteration ", i, ":")
         echo("   Loss = ", formatEng(loss))
         echo("   Predictions = ", A2)

main()
