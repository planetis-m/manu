import times, stats, strformat

type
   Matrix*[T: SomeFloat] = object
      m, n: int # Row and column dimensions.
      data: ptr UncheckedArray[T] # Array for internal storage of elements.

template createData[T](size): ptr UncheckedArray[T] =
   cast[ptr UncheckedArray[T]](alloc(size * sizeof(T)))

template createData0[T](size): ptr UncheckedArray[T] =
   cast[ptr UncheckedArray[T]](alloc0(size * sizeof(T)))

proc `=destroy`*[T](m: var Matrix[T]) =
   if m.data != nil:
      dealloc(m.data)
      m.data = nil
      m.m = 0
      m.n = 0

proc `=`*[T](a: var Matrix[T]; b: Matrix[T]) =
   if a.data != b.data:
      `=destroy`(a)
      a.m = b.m
      a.n = b.n
      if b.data != nil:
         let len = b.m * b.n
         a.data = createData[T](len)
         copyMem(a.data, b.data, len * sizeof(T))

proc matrix*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
   ## Construct an m-by-n matrix of zeros.
   result.m = m
   result.n = n
   result.data = createData0[T](m * n)

proc matrixUninit*[T: SomeFloat](m, n: int): Matrix[T] {.inline.} =
   ## Construct an m-by-n matrix. Note that the matrix will be uninitialized.
   result.m = m
   result.n = n
   result.data = createData[T](m * n)

proc identity*[T: SomeFloat](m: int): Matrix[T] =
   ## Generate identity matrix.
   ##
   ## ``return``: An m-by-m matrix with ones on the diagonal and zeros elsewhere.
   result = matrix[T](m, m)
   for i in 0 ..< m:
      result.data[i * m + i] = T(1.0)

proc identity2*[T: SomeFloat](m: int): Matrix[T] =
   ## Generate identity matrix.
   ##
   ## ``return``: An m-by-m matrix with ones on the diagonal and zeros elsewhere.
   result = matrixUninit[T](m, m)
   for i in 0 ..< m:
      for j in 0 ..< m:
         if i == j:
            result.data[i * m + j] = T(1.0)
         else:
            result.data[i * m + j] = T(0.0)

proc warmup() =
   # Warmup - make sure cpu is on max perf
   let start = cpuTime()
   var a = 123
   for i in 0 ..< 300_000_000:
      a += i * i mod 456
      a = a mod 789
   let dur = cpuTime() - start
   echo &"Warmup: {dur:>4.4f} s"

proc printStats(name: string, stats: RunningStat, dur: float) =
   echo &"""{name}:
   Collected {stats.n} samples in {dur:>4.4f} s
   Average time: {stats.mean * 1000:>4.4f} ms
   Stddev  time: {stats.standardDeviationS * 1000:>4.4f} ms
   Min     time: {stats.min * 1000:>4.4f} ms
   Max     time: {stats.max * 1000:>4.4f} ms"""

template bench(name, samples, code: untyped) =
   proc runBench() {.gensym, nimcall.} =
      var stats: RunningStat
      let globalStart = cpuTime()
      for i in 1 .. samples:
         let start = cpuTime()
         code
         let duration = cpuTime() - start
         stats.push duration
      let globalDuration = cpuTime() - globalStart
      printStats(name, stats, globalDuration)
   runBench()

proc main =
   const m = 1_000
   warmup()
   bench("init0", m):
      let A = identity[float](m)
   bench("no init", m):
      let A = identity2[float](m)

main()
