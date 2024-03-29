# Definitions copy-pasted from Nim/lib/system/memalloc.nim
# It would be nice if these were finally exported.

{.push stackTrace: off.}

include system/bitmasks

template `+!`(p: pointer, s: SomeInteger): pointer =
  cast[pointer](cast[int](p) +% int(s))

template `-!`(p: pointer, s: SomeInteger): pointer =
  cast[pointer](cast[int](p) -% int(s))

proc alignedAlloc(size, align: Natural): pointer =
  if align <= MemAlign:
    when compileOption("threads"):
      result = allocShared(size)
    else:
      result = alloc(size)
  else:
    # allocate (size + align - 1) necessary for alignment,
    # plus 2 bytes to store offset
    when compileOption("threads"):
      let base = allocShared(size + align - 1 + sizeof(uint16))
    else:
      let base = alloc(size + align - 1 + sizeof(uint16))
    # memory layout: padding + offset (2 bytes) + user_data
    # in order to deallocate: read offset at user_data - 2 bytes,
    # then deallocate user_data - offset
    let offset = align - (cast[int](base) and (align - 1))
    cast[ptr uint16](base +! (offset - sizeof(uint16)))[] = uint16(offset)
    result = base +! offset

proc alignedAlloc0(size, align: Natural): pointer =
  if align <= MemAlign:
    when compileOption("threads"):
      result = allocShared0(size)
    else:
      result = alloc0(size)
  else:
    # see comments for alignedAlloc
    when compileOption("threads"):
      let base = allocShared0(size + align - 1 + sizeof(uint16))
    else:
      let base = alloc0(size + align - 1 + sizeof(uint16))
    let offset = align - (cast[int](base) and (align - 1))
    cast[ptr uint16](base +! (offset - sizeof(uint16)))[] = uint16(offset)
    result = base +! offset

proc alignedDealloc(p: pointer, align: int) =
  if align <= MemAlign:
    when compileOption("threads"):
      deallocShared(p)
    else:
      dealloc(p)
  else:
    # read offset at p - 2 bytes, then deallocate (p - offset) pointer
    let offset = cast[ptr uint16](p -! sizeof(uint16))[]
    when compileOption("threads"):
      deallocShared(p -! offset)
    else:
      dealloc(p -! offset)

{.pop.}
