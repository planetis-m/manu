mode = ScriptMode.Verbose

packageName   = "morpheus"
version       = "1.0.1"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "htmldocs", "examples", "experiments"]

requires "nim >= 0.18.0"

switch "forceBuild"

proc configForTests() =
  switch "hints", "off"
  switch "linedir", "on"
  switch "stacktrace", "on"
  switch "linetrace", "on"
  switch "debuginfo"
  switch "path", "."
  switch "run"

task test, "run tests":
  configForTests()
  setCommand "c", "tests/testMatrix.nim"

task docs, "generate documentation":
  exec("mkdir -p htmldocs/morpheus")
  switch "project"
  switch "docSeeSrcUrl", "https://github.com/notTito/morpheus/blob/master"
  setCommand "doc", "morpheus.nim"
