{.passC: "-march=native -ffast-math".}

const
   memAlign = 256

type
   Matrix*[T: SomeFloat] = object
      m, n: int # Row and column dimensions.
      p: ptr Payload[T]

   Payload[T] = object
      pad: byte
      data {.align(memAlign).}: UncheckedArray[T] # Array for internal storage of elements.

proc align(address, alignment: int): int =
   result = (address + (alignment - 1)) and not (alignment - 1)

template createData[T](size): ptr Payload[T] =
   cast[ptr Payload[T]](allocShared(align(sizeof(Payload[T]), memAlign - 16) + size * sizeof(T)))

template createData0[T](size): ptr Payload[T] =
   cast[ptr Payload[T]](allocShared0(align(sizeof(Payload[T]), memAlign - 16) + size * sizeof(T)))

proc `=destroy`*[T](m: var Matrix[T]) =
   if m.p != nil:
      dealloc(m.p)
      m.p = nil
      m.m = 0
      m.n = 0

proc `=`*[T](a: var Matrix[T]; b: Matrix[T]) =
   if a.p != b.p:
      `=destroy`(a)
      a.m = b.m
      a.n = b.n
      if b.p != nil:
         let len = b.m * b.n
         a.p = createData[T](len)
         copyMem(unsafeAddr a.p.data[0], unsafeAddr b.p.data[0], len * sizeof(T))

proc matrix*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
   ## Construct an m-by-n matrix of zeros.
   result.m = m
   result.n = n
   result.p = createData0[T](m * n)

proc matrixUninit*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
   ## Construct an m-by-n matrix. Note that the matrix will be uninitialized.
   result.m = m
   result.n = n
   result.p = createData[T](m * n)

proc `[]`*[T](m: Matrix[T], i, j: int): lent T {.inline.} =
   result = m.p.data[i * m.n + j]

proc `[]`*[T](m: var Matrix[T], i, j: int): var T {.inline.} =
   m.p.data[i * m.n + j]

proc `[]=`*[T](m: var Matrix[T], i, j: int, s: T) {.inline.} =
   m.p.data[i * m.n + j] = s

import strutils

proc main =
   var a = matrix[float](2, 3)
   echo cast[ByteAddress](a.p).toHex
   echo cast[ByteAddress](a.p.data).toHex
   var b = matrix[float](2, 3)
   echo cast[ByteAddress](b.p).toHex
   echo cast[ByteAddress](b.p.data).toHex
   var c = matrix[float](2, 3)
   echo cast[ByteAddress](c.p).toHex
   echo cast[ByteAddress](c.p.data).toHex
   var d = matrix[float](2, 3)
   echo cast[ByteAddress](d.p).toHex
   echo cast[ByteAddress](d.p.data).toHex

main()
