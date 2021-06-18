# See: https://peterroelants.github.io/posts/neural-network-implementation-part04/
# Implementing an Artificial Neural Network in manu
import math, strutils, manu/matrix

proc sigmoid(s: float32): float32 {.inline.} =
  result = 1'f32 / (1'f32 + exp(-s))
makeUniversal(sigmoid)

proc loss(y, t: float32): float32 {.inline.} =
  result = t * ln(y) + (1'f32 - t) * ln(1'f32 - y)
makeUniversalBinary(loss)

proc main =
  const
    m = 4 # batch length
    nodes = 3
    rate = 0.01'f32
  let
    X = matrix(2, @[0'f32, 0, 0, 1, 1, 0, 1, 1])
    Y = matrix(1, @[0'f32, 1, 1, 0])
  var
    # Layer 1
    W1 = randMatrix(2, nodes, -1'f32..1'f32)
    b1 = zeros32(1, nodes)
    # Layer 2
    W2 = randMatrix(nodes, 1, -1'f32..1'f32)
    b2 = 0'f32
  for i in 1 .. 1000:
    let
      # Foward Prop
      # Layer 1
      Z1 = X * W1 + RowVector32(b1) # broadcast bias to (m, nodes)
      A1 = sigmoid(Z1)
      # Layer 2
      Z2 = A1 * W2 + b2 # scalar to (m, 1)
      A2 = sigmoid(Z2)
      # Cross Entropy
      loss = -sum(loss(A2, Y)) / m.float32
      # Back Prop
      # Layer 2
      dZ2 = A2 - Y
      db2 = sum(dZ2)
      dW2 = A1.transpose * dZ2
      # Layer 1
      dZ1 = (dZ2 * W2.transpose) *. (1'f32 - A1) *. A1
      db1 = sumColumns(dZ1)
      dW1 = X.transpose * dZ1
    # Gradient Descent
    # Layer 2
    W2 -= rate * dW2
    b2 -= rate * db2
    # Layer 1
    W1 -= rate * dW1
    b1 -= rate * db1
    if i mod 250 == 0:
      echo(" Iteration ", i, ":")
      echo("   Loss = ", formatEng(loss))
      echo("   Predictions =\n", A2)

main()
