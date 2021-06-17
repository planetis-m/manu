packageName   = "manu"
version       = "2.2.1"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "docs", "examples", "experiments"]

requires "nim >= 1.5.0"
requires "nake"

task test, "run tests":
  exec "nake test"
