# From: https://gist.github.com/Jeraldy/1aa6ae6fefa46b7a9cc02b6573cfeefe
# Implementing an Artificial Neural Network in manu
import math, strutils, manu

proc sigmoid(m: sink Matrix): Matrix =
   result = m
   for i in 0 ..< result.m:
      for j in 0 ..< result.n:
         result[i, j] = 1.0 / (1.0 + exp(-result[i, j]))

proc crossEntropy(batchLen: int, y: Matrix, a: sink Matrix): float =
   var a = a
   for i in 0 ..< a.m:
      for j in 0 ..< a.n:
         a[i, j] = y[i, j] * ln(a[i, j]) + (1.0 - y[i, j]) * ln(1.0 - a[i, j])
   result = -sum(a) / batchLen.float

proc main =
   const
      m = 4
      nodes = 400
   let
      X = matrix(m, @[0.0, 0, 1, 1, 0, 1, 0, 1]) # X.transpose
      Y = matrix(m, @[0.0, 1, 1, 0]) # Y.transpose
   var
      rate = 1.0
      anneal = 0.99
      # LAYER 1
      W1 = randMatrix(nodes, 2)
      b1 = matrix(nodes, m) # zeros
      # LAYER 2
      W2 = randMatrix(1, nodes)
      b2 = matrix(1, m) # zeros
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
         dW2 = (dZ2 * A1.transpose) / m.float
         db2 = dZ2 / m.float
         # LAYER 1
         dZ1 = (W2.transpose * dZ2) .* (1.0 - A1 .* A1)
         dW1 = (dZ1 * X.transpose) / m.float
         db1 = dZ1 / m.float
      # Gradient Descent
      # LAYER 2
      W2 -= rate * dW2
      b2 -= rate * db2
      # LAYER 1
      W1 -= rate * dW1
      b1 -= rate * db1
      # Loss
      let loss = crossEntropy(m, Y, A2)
      rate = rate * anneal
      if i mod 20 == 0:
         echo("Iteration : ", i)
         echo("  Loss = ", loss.formatEng)
         echo("  Predictions = ", A2)

main()
