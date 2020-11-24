# Special case:
# Element-by-element left division, ``C = A.\B``
template `\.`*[T](b: Matrix[T], a: Matrix[T]): Matrix[T] = `/.`(a, b)
template `\.`*[T](b: Matrix[T], a: ColVector[T]): Matrix[T] = `/.`(a, b)
template `\.`*[T](b: Matrix[T], a: RowVector[T]): Matrix[T] = `/.`(a, b)
template `\.`*[T](b: ColVector[T], a: RowVector[T]): Matrix[T] = `/.`(a, b)
template `\.`*[T](b: RowVector[T], a: ColVector[T]): Matrix[T] = `/.`(a, b)
template `\.`*[T](b: ColVector[T], a: Matrix[T]): Matrix[T] = `/.`(a, b)
template `\.`*[T](b: RowVector[T], a: Matrix[T]): Matrix[T] = `/.`(a, b)
# in place
template `\.=`*[T](a: Matrix[T], b: Matrix[T]) = `/.=`(a, b)
template `\.=`*[T](a: Matrix[T], b: ColVector[T]) = `/.=`(a, b)
template `\.=`*[T](a: Matrix[T], b: RowVector[T]) = `/.=`(a, b)
