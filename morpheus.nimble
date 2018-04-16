mode = ScriptMode.Verbose

packageName   = "morpheus"
version       = "0.9.1"
author        = "Antonis Geralis"
description   = "Nim Matrix module"
license       = "MIT"
skipDirs = @["tests", "htmldocs", "examples"]

requires "nim >= 0.18.0"

--forceBuild

proc configForTests() =
  --hints: off
  --linedir: on
  --stacktrace: on
  --linetrace: on
  --debuginfo
  --path: "."
  --run

task test, "run tests":
  configForTests()
  setCommand "c", "tests/testMatrix.nim"

task docs, "generate documentation":
  exec("mkdir -p htmldocs/morpheus")
  --project
  --docSeeSrcUrl: "https://github.com/notTito/morpheus/blob/master"
  setCommand "doc", "morpheus.nim"
